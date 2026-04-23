mod sdust;
mod tr;

use crate::tr::find_repeats;
use clap::{
    CommandFactory, Parser,
    builder::styling::{AnsiColor, Style, Styles},
};

use clio::Input;
use needletail::{
    parse_fastx_reader,
    parser::{FastxReader, SequenceRecord},
};
use std::io::{self, IsTerminal, Write};
use std::ops::RangeInclusive;
use std::process;
use std::str::{Utf8Error, from_utf8};

const FRACTION_RANGE: RangeInclusive<f64> = 0.0..=1.0;

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
    /// (requires --enable-itr-trimming)
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

fn format_record(
    record: &SequenceRecord<'_>,
    sequence: &[u8],
    has_dtr: bool,
    has_itr: bool,
    tr_length: usize,
    include_tr_info: bool,
    disable_trimming: bool,
) -> Result<String, Utf8Error> {
    let header = from_utf8(record.id())?;
    let sequence = from_utf8(sequence)?;
    let trimmed_sequence = if (has_dtr || has_itr) & !disable_trimming {
        &sequence[..sequence.len() - tr_length]
    } else {
        sequence
    };
    let header_line = if include_tr_info {
        match (has_dtr, has_itr) {
            (true, _) => format!(">{} tr=dtr tr_length={}", header, tr_length),
            (_, true) => format!(">{} tr=itr tr_length={}", header, tr_length),
            _ => format!(">{} tr=none tr_length=0", header),
        }
    } else {
        format!(">{}", header)
    };
    let wrapped_sequence = textwrap::fill(trimmed_sequence, 80);
    Ok(format!("{}\n{}", header_line, wrapped_sequence))
}

fn write_record_to_stdout(record: String) {
    match writeln!(io::stdout(), "{}", record) {
        Ok(_) => (),
        Err(e) => match e.kind() {
            io::ErrorKind::BrokenPipe => std::process::exit(0),
            _ => eprintln!("Error writing to stdout: {}", e),
        },
    }
}

fn pipeline(
    mut reader: Box<dyn FastxReader>,
    enable_itr_identification: bool,
    disable_dtr_trimming: bool,
    min_length: usize,
    ignore_low_complexity: bool,
    max_low_complexity_frac: f64,
    ignore_ambiguous: bool,
    max_ambiguous_frac: f64,
    exclude_non_tr_seqs: bool,
    include_tr_info: bool,
    disable_trimming: bool,
) {
    while let Some(record) = reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(e) => {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        };

        let sequence = &record.seq();

        let (has_dtr, has_itr, tr_length) = find_repeats(
            sequence,
            min_length,
            disable_dtr_trimming,
            enable_itr_identification,
            ignore_low_complexity,
            max_low_complexity_frac,
            ignore_ambiguous,
            max_ambiguous_frac,
        );

        if !exclude_non_tr_seqs || has_dtr || has_itr {
            match format_record(
                &record,
                sequence,
                has_dtr,
                has_itr,
                tr_length,
                include_tr_info,
                disable_trimming,
            ) {
                Ok(formatted_record) => write_record_to_stdout(formatted_record),
                Err(e) => eprintln!("Error formatting record: {}", e),
            };
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
                        process::exit(0);
                    }
                    // If stdin is invalid but there are other inputs, skip it
                    continue;
                }
                // If the error is from a file input, report and exit
                eprintln!(
                    "Error: failed to create reader for {}: {}",
                    input_display, error_msg
                );
                process::exit(1);
            }
        };

        pipeline(
            reader,
            cli.enable_itr_identification,
            cli.disable_dtr_trimming,
            cli.min_length,
            cli.ignore_low_complexity,
            cli.max_low_complexity_frac,
            cli.ignore_ambiguous,
            cli.max_ambiguous_frac,
            cli.exclude_non_tr_seqs,
            cli.include_tr_info,
            cli.disable_trimming,
        );
    }
}
