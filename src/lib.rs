//! Functionality for reading and working with recombination maps.
//!
//! [`RecMap`] objects can be created from reading in a HapMap-formatted
//! recombination map. Note that since the HapMap recombination format does
//! not include the chromosome lengths, this must be specified too.
//! A convenience function [`read_seqlens`] is provided to read in TSV-formatted
//! "genome" files of the chromosome names and lengths.
//!
//! Here is a example which loads a recombination map from a HapMap-formatted
//! recombination map and calculates the total map lengths.
//!
//! ```no_run
//! use recmap::prelude::*;
//! let seqlens = read_seqlens("hg38_seqlens.tsv")
//!                   .expect("could not read seqlens");
//! let rec_map = RecMap::from_hapmap("decode_2019_map.txt", &seqlens)
//!                   .expect("cannot read hapmap");
//!
//! for (name, rate_map) in rec_map.iter() {
//!     println!("{}\t{}", name, rate_map.total_map_length().unwrap());
//! }
//! ```
//!
//! This example can be run on the command line with:
//!
//! ```bash
//! cargo run --example  calc_map_lengths --  --seqlens hg38_seqlens.tsv decode_2019_map.txt
//! ```
//! One of the most common tasks when working with recombination maps is to
//! estimate the map position of arbitrary markers, which is usually done by linear
//! interpolation. `RecMap` provides an easy way to do this for one position
//! (`RecMap.interpolate_map_position()`) and for many positions, with
//! `RecMap.interpolate_map_positions()`:
//!
//!
//! ```no_run
//! use recmap::prelude::*;
//! let seqlens = read_seqlens("hg38_seqlens.tsv")
//!                   .expect("could not read seqlens");
//! let rec_map = RecMap::from_hapmap("decode_2019_map.txt", &seqlens)
//!                   .expect("cannot read hapmap");
//!
//! let positions = vec![11975064, 15007450];
//! rec_map.interpolate_map_positions("chr1", &positions);
//!
//! ```

pub mod file;
mod numeric;
pub mod recmap;

pub use recmap::{read_seqlens, RateFloat, RecMap, RecMapError, CM_MB_CONVERSION};

pub mod prelude {
    pub use crate::recmap::{read_seqlens, RateFloat, RecMap, RecMapError, CM_MB_CONVERSION};
}

#[cfg(test)]
mod tests {}
