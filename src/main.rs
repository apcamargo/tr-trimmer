mod sdust;
mod tr;

use crate::tr::{RepeatOptions, RepeatScratch, find_repeats_with_scratch};
use clap::{
    CommandFactory, Parser,
    builder::styling::{AnsiColor, Style, Styles},
};

use clio::Input;
use needletail::{parse_fastx_reader, parser::FastxReader};
use std::io::{self, BufWriter, IsTerminal, Write};
use std::ops::RangeInclusive;
use std::process;
use std::str::{Utf8Error, from_utf8};

const FRACTION_RANGE: RangeInclusive<f64> = 0.0..=1.0;
const FASTA_LINE_WIDTH: usize = 80;

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Cyan.on_default().bold())
    .usage(AnsiColor::Yellow.on_default().bold())
    .literal(AnsiColor::Yellow.on_default().bold())
    .placeholder(Style::new().dimmed());

fn fraction_in_range(s: &str) -> Result<f64, String> {
    let f: f64 = s.parse().map_err(|_| format!("`{s}` isn't a fraction"))?;
    if FRACTION_RANGE.contains(&f) {
        Ok(f)
    } else {
        Err(format!(
            "value should be between {} and {}",
            FRACTION_RANGE.start(),
            FRACTION_RANGE.end()
        ))
    }
}

/// Trim terminal repeats from sequences in FASTA files
#[derive(Parser)]
#[command(version, about, max_term_width = 79, styles = STYLES)]
struct Cli {
    /// Input file(s). Use '-' for stdin
    #[clap(default_value = "-")]
    input: Vec<Input>,

    /// Identify inverted terminal repeats (ITRs) from sequences
    #[clap(
        short = 'i',
        long,
        default_value = "false",
        help_heading = "Terminal repeat identification"
    )]
    enable_itr_identification: bool,

    /// Disable identification of direct terminal repeats (DTRs) from sequences
    /// (requires --enable-itr-identification)
    #[clap(
        short = 'd',
        long,
        requires = "enable_itr_identification",
        default_value = "false",
        help_heading = "Terminal repeat identification"
    )]
    disable_dtr_trimming: bool,

    /// Minimum length of terminal repeat
    #[clap(
        short = 'l',
        long,
        default_value = "21",
        help_heading = "Terminal repeat identification"
    )]
    min_length: usize,

    /// Ignore terminal repeats that contain a high proportion of low complexity
    /// sequences
    #[clap(
        long,
        short = 'c',
        default_value = "false",
        help_heading = "Terminal repeat filtering"
    )]
    ignore_low_complexity: bool,

    /// Maximum fraction of the terminal repeat length that is comprised of
    /// low-complexity sequence
    #[clap(
        long,
        value_parser = fraction_in_range,
        default_value = "0.5",
        requires = "ignore_low_complexity",
        help_heading = "Terminal repeat filtering"
    )]
    max_low_complexity_frac: f64,

    /// Ignore terminal repeats that contain a high proportion of ambiguous
    /// bases (e.g. 'N')
    #[clap(
        long,
        short = 'n',
        default_value = "false",
        help_heading = "Terminal repeat filtering"
    )]
    ignore_ambiguous: bool,

    /// Maximum fraction of the terminal repeat length that is comprised of
    /// ambiguous bases
    #[clap(
        long,
        value_parser = fraction_in_range,
        default_value = "0.0",
        requires = "ignore_ambiguous",
        help_heading = "Terminal repeat filtering"
    )]
    max_ambiguous_frac: f64,

    /// Retain only the sequences for which terminal repeats were identified
    #[clap(short = 'x', long, default_value = "false", help_heading = "Output")]
    exclude_non_tr_seqs: bool,

    /// Add terminal repeat information to the sequence headers (e.g.,
    /// 'tr=dtr tr_length=55')
    #[clap(short = 'a', long, default_value = "false", help_heading = "Output")]
    include_tr_info: bool,

    /// Disable trimming of terminal repeats from sequences. Can be used with
    /// `--include-tr-info` or `--exclude-non-tr-seqs` to identify and report
    /// sequences with terminal repeats without modifying the sequences
    #[clap(short = 't', long, default_value = "false", help_heading = "Output")]
    disable_trimming: bool,
}

fn create_fasta_reader(input: Input) -> Result<Box<dyn FastxReader>, String> {
    if input.can_seek() && input.is_empty() == Some(true) {
        return Err("the input file is empty".to_string());
    }
    parse_fastx_reader(input).map_err(|e| e.to_string())
}

fn write_fasta_record(
    writer: &mut impl Write,
    header: &[u8],
    sequence: &[u8],
    repeat: RepeatResult,
    output_options: FastaOutputOptions,
) -> Result<(), FastaWriteError> {
    from_utf8(header)?;
    let output_length = if repeat.is_terminal() && !output_options.disable_trimming {
        sequence.len() - repeat.length
    } else {
        sequence.len()
    };

    let sequence = from_utf8(sequence)?;
    let wrapped_sequence = textwrap::fill(&sequence[..output_length], FASTA_LINE_WIDTH);

    writer.write_all(b">")?;
    writer.write_all(header)?;
    if output_options.include_tr_info {
        match (repeat.has_dtr, repeat.has_itr) {
            (true, _) => write!(writer, " tr=dtr tr_length={}", repeat.length)?,
            (_, true) => write!(writer, " tr=itr tr_length={}", repeat.length)?,
            _ => writer.write_all(b" tr=none tr_length=0")?,
        }
    }
    writer.write_all(b"\n")?;

    writer.write_all(wrapped_sequence.as_bytes())?;
    writer.write_all(b"\n")?;
    Ok(())
}

#[derive(Clone, Copy, Debug)]
struct RepeatResult {
    has_dtr: bool,
    has_itr: bool,
    length: usize,
}

impl RepeatResult {
    fn is_terminal(self) -> bool {
        self.has_dtr || self.has_itr
    }
}

#[derive(Clone, Copy, Debug)]
struct FastaOutputOptions {
    include_tr_info: bool,
    disable_trimming: bool,
}

#[derive(Debug)]
enum FastaWriteError {
    Io(io::Error),
    Utf8(Utf8Error),
}

impl From<io::Error> for FastaWriteError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<Utf8Error> for FastaWriteError {
    fn from(error: Utf8Error) -> Self {
        Self::Utf8(error)
    }
}

fn flush_writer_or_exit(writer: &mut impl Write) {
    if let Err(error) = writer.flush() {
        handle_output_error(error);
    }
}

fn exit_after_flush(writer: &mut impl Write, code: i32) -> ! {
    flush_writer_or_exit(writer);
    process::exit(code);
}

fn handle_output_error(error: io::Error) -> ! {
    if error.kind() == io::ErrorKind::BrokenPipe {
        process::exit(0);
    }
    eprintln!("Error writing to stdout: {error}");
    process::exit(1);
}

fn pipeline(
    mut reader: Box<dyn FastxReader>,
    writer: &mut impl Write,
    repeat_options: RepeatOptions,
    exclude_non_tr_seqs: bool,
    output_options: FastaOutputOptions,
) {
    let mut repeat_scratch = RepeatScratch::default();

    while let Some(record) = reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(e) => {
                eprintln!("Error: {}", e);
                exit_after_flush(writer, 1);
            }
        };

        let sequence = record.seq();

        let (has_dtr, has_itr, length) =
            find_repeats_with_scratch(sequence.as_ref(), repeat_options, &mut repeat_scratch);
        let repeat = RepeatResult {
            has_dtr,
            has_itr,
            length,
        };

        if exclude_non_tr_seqs && !repeat.is_terminal() {
            continue;
        }

        if let Err(error) = write_fasta_record(
            writer,
            record.id(),
            sequence.as_ref(),
            repeat,
            output_options,
        ) {
            match error {
                FastaWriteError::Io(error) => handle_output_error(error),
                FastaWriteError::Utf8(error) => {
                    eprintln!("Error formatting record: {error}");
                }
            }
        }
    }
}

fn main() {
    let cli = Cli::parse();
    let input_count = cli.input.len();

    // If it's an interactive session with no data piped to stdin and files provided,
    // show help and exit
    if input_count == 1 && cli.input[0].is_std() && io::stdin().is_terminal() {
        Cli::command().print_help().unwrap();
        process::exit(0);
    }

    let stdout = io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    let repeat_options = RepeatOptions {
        min_length: cli.min_length,
        disable_dtr_identification: cli.disable_dtr_trimming,
        enable_itr_identification: cli.enable_itr_identification,
        ignore_low_complexity: cli.ignore_low_complexity,
        max_low_complexity_frac: cli.max_low_complexity_frac,
        ignore_ambiguous: cli.ignore_ambiguous,
        max_ambiguous_frac: cli.max_ambiguous_frac,
    };
    let output_options = FastaOutputOptions {
        include_tr_info: cli.include_tr_info,
        disable_trimming: cli.disable_trimming,
    };

    for input in cli.input {
        let is_std = input.is_std();
        let input_display = input.to_string();

        let reader = match create_fasta_reader(input) {
            Ok(reader) => reader,
            Err(error_msg) => {
                if is_std {
                    // If stdin is invalid and it's the only input, show help and exit
                    if input_count == 1 {
                        Cli::command().print_help().unwrap();
                        exit_after_flush(&mut writer, 0);
                    }
                    // If stdin is invalid but there are other inputs, skip it
                    continue;
                }
                // If the error is from a file input, report and exit
                eprintln!(
                    "Error: failed to create reader for {}: {}",
                    input_display, error_msg
                );
                exit_after_flush(&mut writer, 1);
            }
        };

        pipeline(
            reader,
            &mut writer,
            repeat_options,
            cli.exclude_non_tr_seqs,
            output_options,
        );
    }

    flush_writer_or_exit(&mut writer);
}

#[cfg(test)]
mod tests {
    use super::{FastaOutputOptions, FastaWriteError, RepeatResult, write_fasta_record};

    fn repeat(has_dtr: bool, has_itr: bool, length: usize) -> RepeatResult {
        RepeatResult {
            has_dtr,
            has_itr,
            length,
        }
    }

    fn output_options(include_tr_info: bool, disable_trimming: bool) -> FastaOutputOptions {
        FastaOutputOptions {
            include_tr_info,
            disable_trimming,
        }
    }

    #[test]
    fn writes_wrapped_fasta_record() {
        let mut output = Vec::new();
        let sequence = b"ACGT".repeat(21);
        write_fasta_record(
            &mut output,
            b"seq1",
            &sequence,
            repeat(false, false, 0),
            output_options(false, false),
        )
        .unwrap();

        let expected = format!(">seq1\n{}\n{}\n", "ACGT".repeat(20), "ACGT",);
        assert_eq!(output, expected.into_bytes());
    }

    #[test]
    fn writes_tr_info_and_trims_sequence() {
        let mut output = Vec::new();
        write_fasta_record(
            &mut output,
            b"seq1 description",
            b"ACGTACGT",
            repeat(true, false, 4),
            output_options(true, false),
        )
        .unwrap();

        assert_eq!(
            output,
            b">seq1 description tr=dtr tr_length=4\nACGT\n".to_vec()
        );
    }

    #[test]
    fn skips_non_utf8_records() {
        let mut output = Vec::new();
        let error = write_fasta_record(
            &mut output,
            b"seq\xff",
            b"ACGT\xff",
            repeat(false, false, 0),
            output_options(false, false),
        )
        .unwrap_err();

        assert!(matches!(error, FastaWriteError::Utf8(_)));
        assert!(output.is_empty());
    }

    #[test]
    fn trims_trailing_spaces_like_textwrap_fill() {
        let mut output = Vec::new();
        write_fasta_record(
            &mut output,
            b"seq",
            b"ACGTACGT ",
            repeat(false, false, 0),
            output_options(false, false),
        )
        .unwrap();

        assert_eq!(output, b">seq\nACGTACGT\n".to_vec());
    }
}
