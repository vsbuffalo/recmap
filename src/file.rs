//! Encapsulates plaintext and gzip-compressed file input and output.
//!
//! The [`InputFile`] and [`OutputFile`] abstractions are for working with
//! possible gzip-compressed TSV files.
//!
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::Write;
use std::io::{self, BufWriter};
use std::io::{BufRead, BufReader, Read};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FileError {
    #[error("IO error: {0}")]
    IOError(#[from] io::Error),
}

/// Check if a file is a gzipped by looking for the magic numbers
fn is_gzipped_file(file_path: &str) -> io::Result<bool> {
    let mut file = File::open(file_path)?;
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;

    Ok(buffer == [0x1f, 0x8b])
}

/// Represents an input file.
///
/// This struct is used to handle operations on an input file, such as reading from the file.
/// This abstracts how data is read in, allowing for both plaintext and gzip-compressed input
/// to be read through a common interface.
pub struct InputFile {
    pub filepath: String,
    pub comments: Option<Vec<String>>,
    pub header: Option<String>,
    pub skip_lines: usize,
}

impl InputFile {
    /// Constructs a new `InputFile`.
    ///
    /// # Arguments
    ///
    /// * `filepath` - A string slice that holds the path to the file. If the file extension is
    /// `.gz`, `InputFile` will automatically uncompress the input.
    pub fn new(filepath: &str) -> Self {
        Self {
            filepath: filepath.to_string(),
            comments: None,
            header: None,
            skip_lines: 0,
        }
    }

    /// Opens the file and returns a buffered reader.
    ///
    /// If the file is gzip-compressed (indicated by a ".gz" extension), this method will
    /// automatically handle the decompression.
    ///
    /// # Returns
    ///
    /// A result containing a `BufReader<Box<dyn Read>>` on success, or a `FileError` on failure.
    ///
    pub fn reader(&self) -> Result<BufReader<Box<dyn Read>>, FileError> {
        let file = File::open(self.filepath.clone())?;
        //let is_gzipped_name = self.filepath.ends_with(".gz");
        let is_gzipped = is_gzipped_file(&self.filepath)?;
        let reader: Box<dyn Read> = if is_gzipped {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Ok(BufReader::new(reader))
    }

    /// Collects comment lines and/or a line at the start of the file.
    pub fn collect_metadata(
        &mut self,
        comment: &str,
        header: Option<&str>,
    ) -> Result<bool, FileError> {
        let mut buf_reader = self.reader()?;
        let mut comments = Vec::new();
        let mut line = String::new();

        while buf_reader.read_line(&mut line)? > 0 {
            if line.starts_with(comment) {
                comments.push(line.trim_end().to_string());
                self.skip_lines += 1;
            } else if let Some(header_string) = header {
                if line.starts_with(header_string) {
                    self.header = Some(line.trim_end().to_string());
                    self.skip_lines += 1;
                    // We only handle one header line. If there are more, the
                    // file is *very* poorly formatted. So just let downstream
                    // parsing errors catch this. In the future, we could have a specialized
                    // error.
                    break;
                }
                // break on the first non-header/comment line
                break;
            }
            line.clear();
        }

        self.comments = Some(comments);
        Ok(self.skip_lines > 0)
    }

    /// Method to continue reading after skipping the comment and header lines.
    pub fn continue_reading(&self) -> Result<BufReader<Box<dyn Read>>, FileError> {
        let mut buf_reader = self.reader()?;
        let mut skipped_lines = 0;
        let mut line = String::new();

        // skip the lines that were previously read as comments or header
        while skipped_lines < self.skip_lines {
            buf_reader.read_line(&mut line)?;
            skipped_lines += 1;
            line.clear();
        }
        Ok(buf_reader)
    }
}

/// Represents an output file.
///
/// This struct is used to handle operations on an output file, such as writing to the file.
/// This abstracts writing both plaintext and gzip-compressed files.
pub struct OutputFile {
    pub filepath: String,
    pub header: Option<Vec<String>>,
}

impl OutputFile {
    /// Constructs a new `OutputFile`.
    ///
    /// # Arguments
    ///
    /// * `filepath` - A string slice that holds the path to the file. If the file extension is
    /// `.gz`, `OutputFile` will automatically write gzip-compressed output.
    /// * `header` - An optional vector of strings representing commented header lines to be written to the file.
    pub fn new(filepath: &str, header: Option<Vec<String>>) -> Self {
        Self {
            filepath: filepath.to_string(),
            header,
        }
    }

    /// Opens the file and returns a writer.
    ///
    /// If the file path ends with ".gz", the file is treated as gzip-compressed, and the
    /// function will handle compression automatically. If a header is set, it will be written
    /// to the file.
    ///
    /// # Returns
    ///
    /// A result containing a `Box<dyn Write>` on success, or an `io::Error` on failure.
    pub fn writer(&self) -> Result<Box<dyn Write>, io::Error> {
        let outfile = &self.filepath;
        let is_gzip = outfile.ends_with(".gz");
        let mut writer: Box<dyn Write> = if is_gzip {
            Box::new(BufWriter::new(GzEncoder::new(
                File::create(outfile)?,
                Compression::default(),
            )))
        } else {
            Box::new(BufWriter::new(File::create(outfile)?))
        };
        // write header if one is set
        if let Some(entries) = &self.header {
            for entry in entries {
                writeln!(writer, "#{}", entry)?;
            }
        }
        Ok(writer)
    }
}
