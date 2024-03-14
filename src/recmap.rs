use genomap::{GenomeMap, GenomeMapError};
use indexmap::map::IndexMap;
use ndarray::Array2;
use num_traits::Float;
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::fmt::LowerExp;
use std::io;
use std::io::BufRead;
use std::io::Write;
use thiserror::Error;

use crate::file::OutputFile;
use crate::numeric::recomb_dist_matrix;

use super::file::{FileError, InputFile};
use super::numeric::interp1d;

/// The float type for recombination rates.
pub type RateFloat = f64;

/// The main position type in `recmap`.
///
/// This type is currently an unwrapped [`u32`]. This should handle
/// chromosome lengths for nearly all species. In fact, the only exception
/// known so far is lungfush (*Neoceratodus forsteri*), which has a chromosomes
/// that reaches 5.4Gb (<https://www.nature.com/articles/s41586-021-03198-8l>).
/// The [`u32::MAX`] is 4,294,967,295, i.e. 4.29 Gigabases, which means [`u32`] is
/// just barely suitable for even the largest known chromosome. There is a
/// performance and memory-efficiency tradeoff when using [`u64`] over [`u32`],
/// so [`u32`] is used by default since it handles nearly all cases.
///
/// # Feature support for large chromosomes
///
/// If you are working with data from a species with unusually large chromosomes,
/// you can compile `recmap` using the `--features=big-position` option, which will set
/// the [`Position`] and [`PositionOffset`] to [`u64`] and [`i64`], respectively.
///
/// [`u32::MAX`]: std::u32::MAX
#[cfg(not(feature = "big-position"))]
pub type Position = u32;
#[cfg(feature = "big-position")]
pub type Position = u64;

/// The main *signed* position type in recmap, to represent offsets (e.g.
/// for adjust range coordinates, etc).
#[cfg(not(feature = "big-position"))]
pub type PositionOffset = i32;
#[cfg(feature = "big-position")]
pub type PositionOffset = i64;

pub const CM_MB_CONVERSION: RateFloat = 1e-8;
pub const RATE_PRECISION: usize = 8;
pub const SCI_NOTATION_THRESH: usize = 8;

#[derive(Error, Debug)]
pub enum RecMapError {
    #[error("HapMap parsing error: {0}")]
    HapMapParsingError(#[from] csv::Error),
    #[error("IO error: {0}")]
    IOError(#[from] io::Error),
    #[error("File reading eror: {0}")]
    FileError(#[from] FileError),
    #[error("Missing field")]
    MissingField,
    #[error("Failed to parse a column of a HapMap file")]
    ParseError(String),
    #[error("Improper Rate value, either NaN or negative ({0})")]
    ImproperRate(String),
    #[error("Chromosome key '{0}' does not exist in the recombination map")]
    NoChrom(String),
    #[error("HapMap file not sorted")]
    HapMapNotSorted,
    #[error("Lookup out of bounds ({0}:{1})")]
    LookupOutOfBounds(String, Position),
    #[error("Internal Error")]
    InternalError(String),
    #[error("Recombination map overuns sequence length for {0} ({1} > {2})")]
    LengthMismatch(String, Position, Position),
    #[error("GenomeMap Error: error updating GenomeMap")]
    GenomeMapError(#[from] GenomeMapError),
}

/// Read a tab-delimited *genome file* of sequence (i.e. chromosome) names and their lengths.
pub fn read_seqlens(filepath: &str) -> Result<IndexMap<String, Position>, csv::Error> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(filepath)?;

    let mut seqlens = IndexMap::new();

    #[derive(Debug, Serialize, Deserialize, Default)]
    struct SeqLenEntry {
        chrom: String,
        length: Position,
    }

    for result in rdr.deserialize() {
        let record: SeqLenEntry = result?;
        seqlens.insert(record.chrom, record.length);
    }

    Ok(seqlens)
}

/// Storage and methods for a single chromosome's recombination rates and marker positions.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
pub struct RateMap {
    /// The n+1 genomic positions of the markers, 0-indexed and ending at the sequence length.
    pub ends: Vec<Position>,
    /// The n rates in Morgans, between each genomic position.
    pub rates: Vec<RateFloat>,
    /// The n+1 cumulative map lengths at each genomic position.
    pub map_pos: Vec<RateFloat>,
}

/// An iterator over the elements of a [`RateMap`]
pub struct RateMapIter {
    ends: std::vec::IntoIter<Position>,
    rates: std::vec::IntoIter<RateFloat>,
    map_pos: std::vec::IntoIter<RateFloat>,
}

impl Iterator for RateMapIter {
    type Item = (Position, RateFloat, RateFloat);

    fn next(&mut self) -> Option<Self::Item> {
        match (self.ends.next(), self.rates.next(), self.map_pos.next()) {
            (Some(end), Some(rate), Some(map_pos)) => Some((end, rate, map_pos)),
            _ => None,
        }
    }
}

impl IntoIterator for RateMap {
    type Item = (Position, RateFloat, RateFloat);
    type IntoIter = RateMapIter;

    fn into_iter(self) -> Self::IntoIter {
        RateMapIter {
            ends: self.ends.into_iter(),
            rates: self.rates.into_iter(),
            map_pos: self.map_pos.into_iter(),
        }
    }
}

impl RateMap {
    /// Create a new recombination rate map for a single chromosome.
    pub fn new() -> Self {
        Self {
            ends: Vec::new(),
            rates: Vec::new(),
            map_pos: Vec::new(),
        }
    }

    /// Returns the spans (i.e. widths) in basepairs between each marker.
    pub fn span(&self) -> Vec<Position> {
        self.ends
            .windows(2)
            .map(|pair| {
                assert!(
                    pair[1] >= pair[0],
                    "invalid positions encountered while calculating span: {:?}",
                    &pair
                );
                pair[1] - pair[0]
            })
            .collect()
    }

    /// Returns the cumulative mass between each marker.
    pub fn mass(&self) -> Vec<RateFloat> {
        self.rates
            .iter()
            .zip(self.span().iter())
            .map(|(&rate, &span)| rate * span as RateFloat)
            .collect()
    }

    /// Calculate the cumulative map length in Morgans at each position.
    pub fn calc_cumulative_mass(&mut self) {
        let mass = self.mass();
        let cumulative_sum: Vec<_> = mass
            .iter()
            .scan(0.0, |state, &x| {
                *state += x;
                Some(*state)
            })
            .collect();
        self.map_pos = cumulative_sum;
    }

    /// Calculate the total map length
    pub fn total_map_length(&self) -> Option<&RateFloat> {
        self.map_pos.last()
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct RecMap {
    pub map: GenomeMap<RateMap>,
    pub metadata: Option<Vec<String>>,
}

impl RecMap {
    /// Create a new [`RecMap`] from a (possibly gzip-compressed) HapMap-formatted recombination map file.
    ///
    /// The HapMap recombination map format is (rather unfortunately) very poorly
    /// specified. This parser is quite permissive, and skips through comment lines
    /// that begin with `#` and a possible header line that beings with `Chr`.
    /// Note that this parser *does not* read the cumulative map positions. Instead,
    /// the cumulative map positions can be calculated directly from the rates and
    /// the distances between markers.
    ///
    /// The HapMap recombination format looks like:
    ///
    /// ```text
    /// Chromosome      Position(bp)    Rate(cM/Mb)     Map(cM)
    /// chr1    55550   2.981822        0.000000
    /// chr1    82571   2.082414        0.080572
    /// chr1    88169   2.081358        0.092229
    /// chr1    254996  3.354927        0.439456
    /// chr1    564598  2.887498        1.478148
    /// chr1    564621  2.885864        1.478214
    /// chr1    565433  2.883892        1.480558
    /// chr1    568322  2.887570        1.488889
    /// chr1    568527  2.895420        1.489481
    /// ```
    ///
    /// The HapMap format is well-described in the [tskit documentation](https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap).
    ///
    /// # Warnings
    ///
    /// Given that the HapMap recombination map format is poorly specified, and users often
    /// implement their own versions, it is **highly recommended** that you always validate
    /// the parsed recombination map. Ideally check that:
    ///
    ///  - The total map length (i.e. in Morgans or centiMorgans) makes sense. Off-by-one
    ///    errors due to the format not following the 0-indexed end-exclusive format assumed
    ///    by [`RecMap.from_hapmap`] will often lead to obviously erroneous total map lengths.
    ///  - Visually plot the recombination map, looking for outlier rates.
    ///  - If the recombination map includes a fourth
    ///    
    ///
    pub fn from_hapmap(
        filepath: &str,
        seqlens: &IndexMap<String, Position>,
    ) -> Result<RecMap, RecMapError> {
        let mut input_file = InputFile::new(filepath);

        let _has_metadata = input_file.collect_metadata("#", Some("Chr"))?;
        let reader = input_file.continue_reading()?;

        let mut rec_map: GenomeMap<RateMap> = GenomeMap::new();
        let mut last_chrom: Option<String> = None;
        let mut last_end: Option<Position> = None;

        // this is used for validation only
        let mut has_fourth_column = false;
        let mut map_positions: Vec<RateFloat> = Vec::new();

        for result in reader.lines() {
            let line = result?;
            let fields: Vec<&str> = line.split_whitespace().collect();

            // get the chrom column
            let chrom = fields.first().ok_or(RecMapError::MissingField)?.to_string();

            // Our parser will always add on the chromosome end. Most times the
            // HapMap file will *not* have this, but we still check.
            match last_chrom {
                None => {}
                Some(last) => {
                    if chrom != last {
                        // We're on a new chromosome. Insert the last entry into rec_map
                        if let Some(seq_len) = seqlens.get(&last) {
                            let chrom_entry = rec_map.entry_or_default(&last);

                            // Assuming the rate is 0 for these extra entries
                            if let Some(&last_end) = chrom_entry.ends.last() {
                                if last_end != *seq_len {
                                    chrom_entry.ends.push(*seq_len);
                                }
                            }
                        }
                        // reset the end too
                        last_end = None;
                    }
                }
            }

            // Update last_chrom and last_end
            last_chrom = Some(chrom.clone());

            // get the position and rate column, parsing into proper numeric types
            let end_str = fields.get(1).ok_or(RecMapError::MissingField)?;
            let end: Position = end_str.parse().map_err(|_| {
                RecMapError::ParseError(format!("Failed to parse end from string: {}", end_str))
            })?;

            // Check that everything is sorted
            match last_end {
                None => Some(end),
                Some(last_end) => {
                    if end >= last_end {
                        return Err(RecMapError::HapMapNotSorted);
                    }
                    Some(end)
                }
            };

            let rate_str = fields.get(2).ok_or(RecMapError::MissingField)?;
            let rate: RateFloat = rate_str.parse().map_err(|_| {
                RecMapError::ParseError(format!("Failed to parse rate from string: {}", rate_str))
            })?;

            // check rate isn't NaN or negative
            if rate.is_nan() || rate < 0.0 {
                return Err(RecMapError::ImproperRate(format!("{}:{}", chrom, end)));
            }

            // if there is a fourth column (total map length) parse it
            if let Some(map_pos_str) = fields.get(3) {
                has_fourth_column = true;
                let map_pos: RateFloat = map_pos_str.parse().map_err(|_| {
                    RecMapError::ParseError(format!(
                        "Failed to parse map position from string: {}",
                        map_pos_str
                    ))
                })?;
                map_positions.push(map_pos);
            }

            // HapMap rates are *always* in cM/Mb, but the natural unit is Morgans, so
            // we convert here.
            let rate = CM_MB_CONVERSION * rate;

            // Insert into GenomeMap, making a new one for this chromosome if needed.
            // Note that we will also pad the first entry so that it's position zero, with rate
            // zero.
            if let Some(chrom_entry) = rec_map.get_mut(&chrom) {
                chrom_entry.ends.push(end);
                chrom_entry.rates.push(rate);
            } else {
                let mut new_rate_map = RateMap::new();

                // If doesn't start with zero, add zero.
                if end != 0 {
                    new_rate_map.ends.push(0);
                    new_rate_map.rates.push(0.0);
                }

                // then add the current entry too
                new_rate_map.ends.push(end);
                new_rate_map.rates.push(rate);

                // if there is a fourth column, we could use it for validation
                // TODO
                if has_fourth_column {
                    //new_rate_map.calc_cumulative_mass();
                    //assert_floats_eq(&new_rate_map.map_pos, map_positions.as_slice(), 0.01);
                    //map_positions.clear();
                }

                rec_map.insert(&chrom, new_rate_map)?;
            }
        }

        // Insert the final entry for the last chromosome outside the loop if needed
        if let Some(last) = last_chrom {
            if let Some(seq_len) = seqlens.get(&last) {
                let chrom_entry = rec_map.get_mut(&last).expect(
                    "internal error: please report at http://github.com/vsbuffalo/recmap/issues",
                );

                // Assuming the rate is 0 for these extra entries
                if chrom_entry.ends.is_empty() || chrom_entry.ends.last().unwrap() != seq_len {
                    if chrom_entry.ends.last().unwrap() >= seq_len {
                        let last_end = *chrom_entry.ends.last().unwrap();
                        return Err(RecMapError::LengthMismatch(last, last_end, *seq_len));
                    }
                    chrom_entry.ends.push(*seq_len);
                }
            }
        }

        let metadata = input_file.comments;
        let mut rec_map = RecMap {
            map: rec_map,
            metadata,
        };
        // generate the map positions from the marker positions
        // and the per-basepair rates..
        rec_map.generate_map_positions();
        Ok(rec_map)
    }

    /// Return the number of chromosomes in the recombination map.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Return if the recombination map is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate over chromosome name and [`RateMap`] tuples.
    pub fn iter(&self) -> impl Iterator<Item = (&String, &RateMap)> {
        self.map.iter()
    }

    /// Generate the cumulative map positions via the marginal recombination rates.
    fn generate_map_positions(&mut self) {
        for rate_map in self.map.values_mut() {
            rate_map.calc_cumulative_mass();
        }
    }

    /// Interpolate the recombination map position at the specified physical position.
    /// This uses linear interpolation.
    ///
    /// # Arguments
    ///  * `name`: the chromosome name.
    ///  * `position`: the physical position to estimate the recombination map position at.
    pub fn interpolate_map_position(
        &self,
        name: &str,
        position: Position,
    ) -> Result<RateFloat, RecMapError> {
        let rate_map = self
            .map
            .get(name)
            .ok_or(RecMapError::NoChrom(name.to_string()))?;
        let ends = &rate_map.ends;
        let interp_result = interp1d(&ends[0..ends.len() - 1], 
                                     &rate_map.map_pos, 
                                     position);
        let interpolated_map_pos =
            interp_result.ok_or(RecMapError::LookupOutOfBounds(name.to_string(), position))?;
        Ok(interpolated_map_pos)
    }

    /// Interpolate the recombination map position at the specified physical positions.
    /// This uses linear interpolation.
    ///
    /// # Arguments
    ///  * `name`: the chromosome name.
    ///  * `position`: the physical position to estimate the recombination map position at.
    pub fn interpolate_map_positions(
        &self,
        chrom: &str,
        positions: &[Position],
    ) -> Result<Vec<RateFloat>, RecMapError> {
        let positions: Vec<RateFloat> = positions
            .iter()
            .map(|p| self.interpolate_map_position(chrom, *p))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(positions)
    }

    /// Build the pairwise recombination distance matrix for the specified chromosome.
    ///
    /// Creates a `positions_x.len() x positions_y.len()` matrix of recombination
    /// *distances* (in Morgans), for the supplied set of positions on the physical
    /// map.
    ///
    /// # Arguments
    ///  * `positions_x`: the first set of marker positions.
    ///  * `positions_y`: the second set of marker positions (just repeat `positions_x` for a
    ///    symmetric distance matrix).
    ///  * `haldane`: whether to convert the recombination distances in *Morgans* to a
    ///      unit-less recombination *fraction*.
    ///  * `rec_floor`: an optional *floor* value; all elements in the matrix less than
    ///      this value will be set to this value. This is sometimes useful in downstream
    ///      processing when zero values create problems.
    ///
    pub fn recomb_dist_matrix(
        &self,
        chrom: &str,
        positions_x: &[Position],
        positions_y: &[Position],
        haldane: bool,
        rec_floor: Option<RateFloat>,
    ) -> Result<Array2<RateFloat>, RecMapError> {
        let x_pos = self.interpolate_map_positions(chrom, positions_x)?;
        let y_pos = self.interpolate_map_positions(chrom, positions_y)?;
        Ok(recomb_dist_matrix(&x_pos, &y_pos, haldane, rec_floor))
    }

    /// Write recombination map to HapMap-formatted file.
    ///
    /// This file has the usual HapMap recombination map header, and columns:
    ///  1. Chromosome name
    ///  2. Position (0-indexed and right-exclusive)
    ///  3. Rate
    ///  4. Map position
    ///
    /// # Arguments
    ///  * `filepath`: The filepath to write the recombination map to. If the filepath
    ///  has an `.gz` extension, the output will be gzip compressed.
    ///  If `filepath` is `None`, uncompressed output will be written to standard out.
    pub fn write_hapmap(&self, filepath: Option<&str>) -> Result<(), RecMapError> {
        let mut writer: Box<dyn Write> = match filepath {
            Some(path) => {
                let file = OutputFile::new(path, None);
                file.writer()?
            }
            None => {
                // Use stdout if filepath is None
                Box::new(std::io::stdout())
            }
        };

        // write that weird HapMap header
        writeln!(writer, "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)")?;

        for (chrom, rate_map) in self.map.iter() {
            // write the rows
            let n = rate_map.ends.len();
            // TODO cut off end
            let ends = &rate_map.ends[1..n - 1];

            for (i, end) in ends.iter().enumerate() {
                let rate = rate_map.rates[i + 1];
                let map_pos = rate_map.map_pos[i + 1];

                // Write the record to the file
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    chrom,
                    end,
                    format_float(rate / CM_MB_CONVERSION),
                    format_float(map_pos),
                )?;
            }
        }

        Ok(())
    }

    /// Write recombination map to a BED-like TSV file.
    ///
    /// This file has columns:
    ///  1. Chromosome name
    ///  2. Start position (0-indexed)
    ///  3. End position (0-indexed and right-exclusive)
    ///  4. Rate
    ///
    /// # Arguments
    ///  * `filepath`: The filepath to write the recombination map to. If the filepath
    ///  has an `.gz` extension, the output will be gzip compressed.
    ///  If `filepath` is `None`, uncompressed output will be written to standard out.
    pub fn write_tsv(&self, filepath: Option<&str>) -> Result<(), RecMapError> {
        let mut writer: Box<dyn Write> = match filepath {
            Some(path) => {
                let file = OutputFile::new(path, None);
                file.writer()?
            }
            None => {
                // Use stdout if filepath is None
                Box::new(std::io::stdout())
            }
        };

        for (chrom, rate_map) in self.map.iter() {
            // get the (start, end) ranges from the end points.
            let ranges: Vec<(Position, Position)> = rate_map
                .ends
                .windows(2)
                .map(|pair| (pair[0], pair[1]))
                .collect();

            // write the rows
            for (i, range) in ranges.iter().enumerate() {
                let rate = rate_map.rates[i];

                // Write the record to the file
                let rate_rescaled: RateFloat = rate / CM_MB_CONVERSION;

                let formatted_rate = format!("{:.1$e}", rate_rescaled, RATE_PRECISION);

                // Write the record to the file
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    chrom, range.0, range.1, formatted_rate
                )?;
            }
        }

        Ok(())
    }
}

pub fn format_float<T>(x: T) -> String
where
    T: Float + LowerExp + Display,
{
    let min = T::from(SCI_NOTATION_THRESH).unwrap();
    let max = T::from(SCI_NOTATION_THRESH).unwrap();
    if x.abs().log10() < -min || x.abs().log10() > max {
        format!("{:.1$e}", x, RATE_PRECISION)
    } else {
        format!("{:.*}", RATE_PRECISION, x)
    }
}

#[cfg(test)]
mod tests {
    use super::Position;
    use crate::{
        numeric::{assert_float_eq, assert_floats_eq},
        prelude::*,
    };
    use indexmap::IndexMap;
    use tempfile::tempdir;

    fn mock_seqlens() -> IndexMap<String, Position> {
        let seqlens = indexmap::indexmap! {
            "chr1".to_string() => 25,
            "chr2".to_string() => 32,
            "chr3".to_string() => 22,
        };
        seqlens
    }

    fn read_hapmap() -> RecMap {
        let seqlens = mock_seqlens();
        let rec_map = RecMap::from_hapmap("tests_data/test_hapmap.txt", &seqlens).unwrap();
        rec_map
    }

    fn to_morgans(x: Vec<RateFloat>) -> Vec<RateFloat> {
        x.iter().map(|v| v * CM_MB_CONVERSION).collect()
    }

    #[test]
    fn test_read_hapmap() {
        let rm = read_hapmap();
        assert_eq!(rm.len(), 3);
        assert!(!rm.is_empty());

        dbg!(&rm.map.get("chr1").unwrap().map_pos);

        let chr1_map = rm.map.get("chr1").unwrap();
        assert_eq!(chr1_map.ends.len(), chr1_map.rates.len() + 1);
        assert_eq!(chr1_map.ends, vec![0, 10, 15, 20, 25]);
        assert_eq!(chr1_map.rates, to_morgans(vec![0.0, 1.10, 1.50, 4.33]));

        // cumulative map calculation:
        // [0-10): *implied 0.0*
        // [10-15): 1.10
        // [15-20): 1.50
        // [20-25): 4.33 *to end*
        // ----
        // 5bp * 1.10 = 5.5
        // 5bp * 1.5 = 7.5 + 5.5 = 13
        // 5bp * 4.33 = 21.65 + 13 = 34.65
        // total = 34.65

        let total_len = *chr1_map.total_map_length().unwrap();
        assert_float_eq(total_len, 34.65 * CM_MB_CONVERSION, 1e-3);
        assert_floats_eq(
            &chr1_map.map_pos,
            &to_morgans(vec![0.0, 5.5, 13., 34.65]),
            1e-3,
            );
    }

    #[test]
    #[ignore = "for debugging"]
    fn test_write_hapmap_local() {
        // this writes the output "locally" in the project directory
        // for easier debugging.
        let seqlens = mock_seqlens();

        let rm = read_hapmap();

        let filepath = "test_hapmap.txt";

        rm.write_hapmap(Some(filepath)).unwrap();

        let rm_readin = RecMap::from_hapmap(filepath, &seqlens).unwrap();

        assert_eq!(rm_readin, rm);
    }

    #[test]
    fn test_write_hapmap() {
        let seqlens = mock_seqlens();

        let rm = read_hapmap();

        // temp dir
        let dir = tempdir().unwrap();
        let binding = dir.path().join("test_hapmap.txt");
        let filepath = binding.to_str().unwrap();

        rm.write_hapmap(Some(filepath)).unwrap();

        let rm_readin = RecMap::from_hapmap(filepath, &seqlens).unwrap();

        assert_eq!(rm_readin, rm);
    }
}
