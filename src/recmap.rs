use csv::ReaderBuilder;
use genomap::{GenomeMap, GenomeMapError};
use indexmap::map::IndexMap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::io;
use std::io::Write;
use thiserror::Error;

use crate::file::OutputFile;

use super::file::{FileError, InputFile};
use super::numeric::interp1d;

/// The float type for recombination rates.
pub type RateFloat = f32;

/// The integer type for genomic positions.
///
/// # Developer Notes
/// In the future, applications needing support for chromosomes longer than
/// `i32::MAX` could change this through a `--feature`.
pub type Position = u64;

const CM_MB_CONVERSION: RateFloat = 1e-8;

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
    #[error("Chromosome key '{0}' does not exist")]
    NoChrom(String),
    #[error("HapMap file not sorted")]
    HapMapNotSorted,
    #[error("Lookup out of bounds ({0}:{1})")]
    LookupOutOfBounds(String, Position),
    #[error("Internal Error")]
    InternalError(String),
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
#[derive(Debug, Serialize, Deserialize, Default)]
pub struct RateMap {
    /// The n+1 genomic positions of the markers, 0-indexed and ending at the sequence length.
    pub ends: Vec<Position>,
    /// The n rates in Morgans, between each genomic position.
    pub rates: Vec<RateFloat>,
    /// The n+1 cumulative map lengths at each genomic position.
    pub map_pos: Vec<RateFloat>,
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
        self.ends.windows(2).map(|pair| pair[1] - pair[0]).collect()
    }

    /// Returns the cumulative mass between each marker.
    pub fn mass(&self) -> Vec<RateFloat> {
        self.rates
            .iter()
            .zip(self.span().iter())
            .map(|(&x, &y)| x * y as RateFloat)
            .collect()
    }

    /// Calculate the cumulative map length at each position.
    pub fn calc_cumulative_mass(&mut self) {
        let mass = self.mass();
        //println!("mass: {:?}", mass);
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

pub struct RecMap {
    pub map: GenomeMap<RateMap>,
}

impl RecMap {
    /// Create a new [`RecMap`] from a HapMap-formatted recombination map file.
    ///
    /// This method also supports reading directly from a gzip-compressed file.
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
    /// This parser is a bit more permissive -- it will allow no headers.
    ///
    pub fn from_hapmap(
        filepath: &str,
        seqlens: IndexMap<String, Position>,
    ) -> Result<RecMap, RecMapError> {
        let input_file = InputFile::new(filepath);

        // read one line to check for headers
        let has_header = input_file.has_header("Chromosome")?;

        let buf_reader = input_file.reader()?;

        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(has_header)
            .from_reader(buf_reader);

        let mut rec_map: GenomeMap<RateMap> = GenomeMap::new();

        let mut last_chrom: Option<String> = None;
        let mut last_end: Option<Position> = None;

        for result in rdr.records() {
            let record = result.map_err(RecMapError::HapMapParsingError)?;

            // remove comment lines (TODO could store in metadata)
            if record.get(0).map_or(false, |s| s.starts_with('#')) {
                continue;
            }

            // get the chrom column
            let chrom = record.get(0).ok_or(RecMapError::MissingField)?.to_string();

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
                                    chrom_entry.rates.push(0.0);
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
            let end_str = record.get(1).ok_or(RecMapError::MissingField)?;
            let end: Position = end_str.parse().map_err(|_| {
                dbg!(&end_str);
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

            let rate_str = record.get(2).ok_or(RecMapError::MissingField)?;
            let rate: RateFloat = rate_str.parse().map_err(|_| {
                RecMapError::ParseError(format!("Failed to parse rate from string: {}", rate_str))
            })?;

            // check rate isn't NaN or negative
            if rate.is_nan() || rate < 0.0 {
                return Err(RecMapError::ImproperRate(format!("{}:{}", chrom, end)));
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

                rec_map.insert(&chrom, new_rate_map)?;
            }
        }

        // Insert the final entry for the last chromosome outside the loop if needed
        if let Some(last) = last_chrom {
            if let Some(seq_len) = seqlens.get(&last) {
                let chrom_entry = rec_map.entry_or_default(&last);

                // Assuming the rate is 0 for these extra entries
                if chrom_entry.ends.is_empty() || chrom_entry.ends.last().unwrap() != seq_len {
                    chrom_entry.ends.push(*seq_len);
                }
            }
        }

        let mut rec_map = RecMap { map: rec_map };
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
        let interp_result = interp1d(&ends[0..ends.len() - 1], &rate_map.map_pos, position);
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
    ) -> Result<Array1<RateFloat>, RecMapError> {
        let positions: Vec<RateFloat> = positions
            .iter()
            .map(|p| self.interpolate_map_position(chrom, *p))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Array1::from_vec(positions))
    }

    /// Write recombination map to a BED-like TSV file.
    ///
    /// # Arguments
    ///  * `filepath`: The filepath to write the recombination map to. If the filepath
    ///  has an `.gz` extension, the output will be gzip compressed.
    ///  If `filepath` is `None`, uncompressed output will be written to standard out.
    pub fn write_tsv(&self, filepath: Option<&str>) -> Result<(), RecMapError> {
        let precision = 8;
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

                let formatted_rate = format!("{:.1$e}", rate_rescaled, precision);

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

#[cfg(test)]
mod tests {
    use crate::RecMap;
    use tempfile::tempdir;

    fn read_hapmap() -> RecMap {
        let seqlens = indexmap::indexmap! {
            "chr1".to_string() => 6980669,
            "chr2".to_string() => 6004443,
            "chr3".to_string() => 6026894,
        };
        let rec_map = RecMap::from_hapmap("tests/data/decode_2010_test_map.txt", seqlens).unwrap();

        let dir = tempdir().unwrap();
        let output_path = dir.path().join("output.tsv");

        rec_map
            .write_tsv(Some(output_path.to_str().unwrap()))
            .unwrap();

        dir.close().unwrap();
        rec_map
    }

    #[test]
    fn test_hapmap_read() {
        let rm = read_hapmap();
        assert_eq!(rm.len(), 3);
        assert!(!rm.is_empty());
    }
}
