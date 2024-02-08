use clap::Parser;
use recmap::prelude::*;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Path to the sequence lengths file
    #[clap(long, value_parser)]
    seqlens: String,

    /// Path to the HapMap recombination map file
    #[clap(value_parser)]
    hapmap: String,
}

fn main() -> Result<(), RecMapError> {
    let args = Args::parse();
    let seqlens = read_seqlens(&args.seqlens)?;
    let rec_map = RecMap::from_hapmap(&args.hapmap, seqlens)?;

    for (name, rate_map) in rec_map.iter() {
        println!("{}\t{}", name, rate_map.total_map_length().unwrap());
    }

    Ok(())
}
