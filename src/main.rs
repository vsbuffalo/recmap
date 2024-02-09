use clap::{Parser, Subcommand};
use recmap::recmap::format_float;
use recmap::recmap::Position;
use recmap::RateFloat;
use recmap::CM_MB_CONVERSION;
use recmap::{read_seqlens, RecMap, RecMapError};
use std::io;
use std::io::BufRead;
use std::io::Write;

const INFO: &str = "\
recmap: manipulate recombination maps
usage: recmap [--help] <subcommand>

Subcommands:
  
  interp: linearly interpolate the positions in a BED file.
 
";

#[derive(Parser)]
#[clap(name = "sdf")]
#[clap(about = INFO)]
struct Cli {
    #[arg(short, long, action = clap::ArgAction::Count)]
    debug: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Linearly interpolate the recombination map position for the specified BED3 file.
    ///
    /// This will output a TSV with the following columns:
    ///
    ///  - chromosome name
    ///  - start position                     (0-indexed, inclusive)
    ///  - end position                       (0-indexed, exclusive)
    ///  - cumulative map position at start   (in centiMorgans)
    ///  - cumulative map position at end     (in centiMorgans)
    ///  - recombination rate in [start, end) (in cM/Mb)
    ///
    /// Example:
    ///
    ///  $ recmap interp --seqlens hg38_seqlens.tsv --hapmap decode_2019_map.txt \
    ///      hg38_1Mb_windows.bed --output decode_2019_map_1Mb_summaries.tsv --header
    ///
    /// It is highly advised all recombination map output is validated visually,
    /// since it is not uncommon for recombination maps to implement non-standard
    /// variants of the HapMap recombination map format.
    Interp {
        /// a TSV file of chromosome names and their lengths
        #[arg(long, required = true)]
        seqlens: String,
        /// the input recombination map in HapMap format
        #[arg(long, required = true)]
        hapmap: String,
        /// the output file path (if not set, uses standard out)
        #[arg(long)]
        output: Option<String>,
        /// a BED-like TSV file of chromosome, position columns to interpolate at
        #[arg(required = true)]
        bedfile: String,
        /// Include a header
        #[arg(long, default_value_t = false)]
        header: bool,
    },
}

fn write_entry(
    writer: &mut Box<dyn Write>,
    recmap: &RecMap,
    chrom: &str,
    starts: &[Position],
    ends: &[Position],
) -> Result<(), RecMapError> {
    let map_starts = recmap.interpolate_map_positions(chrom, starts)?;
    let map_ends = recmap.interpolate_map_positions(chrom, ends)?;

    let iter = starts.iter().zip(ends).zip(map_starts).zip(map_ends);

    for (((start, end), map_start), map_end) in iter {
        let width = end - start;
        let diff = (map_end - map_start) / width as RateFloat;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom,
            start,
            end,
            format_float(map_start * 100.),
            format_float(map_end * 100.),
            format_float(diff / CM_MB_CONVERSION)
        )?;
    }
    Ok(())
}

fn interpolate_map(
    bedfile: &str,
    hapmap: &str,
    seqlens: &str,
    output: Option<&str>,
    header: bool,
) -> Result<(), RecMapError> {
    let sl = read_seqlens(seqlens)?;
    let rm = RecMap::from_hapmap(hapmap, &sl)?;

    // open writer, possibly to stdout
    let mut writer = if let Some(filepath) = output {
        let output = recmap::file::OutputFile::new(filepath, None);
        output.writer()?
    } else {
        Box::new(io::stdout())
    };

    if header {
        // write the header
        writeln!(writer, "chrom\tstart\tend\tstart_map\tend_map\trate")?;
    }

    // open reader
    let mut tsvfile = recmap::file::InputFile::new(bedfile);
    // remove the header, if there
    let _has_metadata = tsvfile.collect_metadata("#", None);
    let reader = tsvfile.continue_reading()?;

    // process input
    let mut current_chrom = None;
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    for result in reader.lines() {
        let line = result?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        let chrom = fields.first().ok_or(RecMapError::MissingField)?.to_string();
        let start_str = fields.get(1).ok_or(RecMapError::MissingField)?;
        let start: Position = start_str.parse().map_err(|_| {
            RecMapError::ParseError(format!("Failed to parse end from string: {}", start_str))
        })?;
        let end_str = fields.get(2).ok_or(RecMapError::MissingField)?;
        let end: Position = end_str.parse().map_err(|_| {
            RecMapError::ParseError(format!("Failed to parse end from string: {}", end_str))
        })?;

        match current_chrom {
            None => current_chrom = Some(chrom),
            Some(ref cur_chrom) => {
                if *cur_chrom != chrom {
                    // we've hit a new chromosome. Run the interpolation.
                    write_entry(&mut writer, &rm, &cur_chrom, &starts, &ends)?;
                    starts.clear();
                    ends.clear();
                    current_chrom = Some(chrom);
                }
                starts.push(start);
                ends.push(end);
            }
        }
    }

    // process last chromosome
    if let Some(cur_chrom) = current_chrom {
        write_entry(&mut writer, &rm, &cur_chrom, &starts, &ends)?;
    }
    Ok(())
}

fn run() -> Result<(), RecMapError> {
    let cli = Cli::parse();
    match &cli.command {
        Some(Commands::Interp {
            bedfile,
            hapmap,
            seqlens,
            output,
            header,
        }) => interpolate_map(bedfile, hapmap, seqlens, output.as_deref(), *header),
        None => {
            println!("{}\n", INFO);
            std::process::exit(1);
        }
    }
}

fn main() {
    match run() {
        Ok(_) => {}
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}
