![Crates.io](https://img.shields.io/crates/v/recmap) ![Crates.io](https://img.shields.io/crates/d/recmap) [![docs](https://docs.rs/recmap/badge.svg)](https://docs.rs/recmap) ![Rust CI](https://github.com/vsbuffalo/recmap/actions/workflows/rust.yml/badge.svg)


# RecMap library (and command line tool) for reading and working with recombination maps in Rust

A `RecMap` object can be created from reading in a HapMap-formatted 
recombination map. Note that since the HapMap recombination format does
not include the chromosome lengths, this must be specified too.
A convenience function `read_seqlens` is provided to read in TSV-formatted
"genome" files of the chromosome names and lengths.

Here is a example which loads a recombination map from a HapMap-formatted 
recombination map and calculates the total map lengths.

```rust
use recmap::prelude::*;
let seqlens = read_seqlens("hg38_seqlens.tsv")
                  .expect("could not read seqlens");
let rec_map = RecMap::from_hapmap("decode_2019_map.txt", seqlens)
                  .expect("cannot read hapmap");

for (name, rate_map) in rec_map.iter() {
    println!("{}\t{}", name, rate_map.total_map_length().unwrap());
}
```

This example can be run on the command line with:

```bash
cargo run --example  calc_map_lengths --  --seqlens hg38_seqlens.tsv decode_2019_map.txt
```

One of the most common tasks when working with recombination maps is to
estimate the map position of arbitrary markers, which is usually done by linear
interpolation. `RecMap` provides an easy way to do this for one position
(`RecMap.interpolate_map_position()`) and for many positions, with 
`RecMap.interpolate_map_positions()`:

```rust
use recmap::prelude::*;
let seqlens = read_seqlens("hg38_seqlens.tsv")
                  .expect("could not read seqlens");
let rec_map = RecMap::from_hapmap("decode_2019_map.txt", seqlens)
                  .expect("cannot read hapmap");

let positions = vec![11975064, 15007450];
rec_map.interpolate_map_positions("chr1", &positions);

```

## Command line tool

Additionally, `recmap` had an optional command line tool feature that
interpolates recombination map positions and recombination rates, given BED3
input:

```
$ recmap interp --seqlens hg38_seqlens.tsv --hapmap decode_2019_map.txt \
   hg38_1Mb_windows.bed --output decode_2019_map_1Mb_summaries.tsv --header
```

Currently the command line tool only has one subcommand, though more features
may be added. Please [file an
issue](https://github.com/vsbuffalo/recmap/issues) if there is a feature you'd
like!

## Installation 

To use the library in your own Rust projects, install with:

```
$ cargo add recmap
```

To install the command line tool, use:

```
$ cargo install recmap --features=cli
```

