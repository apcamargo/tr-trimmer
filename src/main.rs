mod sdust;
mod tr;
use std::io::{self, Write};

use crate::tr::find_repeats;
use clap::Parser;
use clio::Input;
use needletail::parser::SequenceRecord;
use needletail::{parse_fastx_file, parse_fastx_stdin};
use std::ops::RangeInclusive;
use std::process;
use std::str::{from_utf8, Utf8Error};

const FRACTION_RANGE: RangeInclusive<f64> = 0.0..=1.0;

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
#[command(version, about, max_term_width = 79)]
struct Cli {
    /// Input file(s). Use '-' for stdin
    #[clap(value_parser, default_value = "-")]
    input: Vec<Input>,

    /// Identify inverted terminal repeats (ITRs) from sequences
    #[clap(
        short = 'i',
        long,
        value_parser,
        default_value = "false",
        help_heading = "Terminal repeat identification"
    )]
    enable_itr_identification: bool,

    /// Disable identification of direct terminal repeats (DTRs) from sequences
    /// (requires --enable-itr-trimming)
    #[clap(
        short = 'd',
        long,
        value_parser,
        requires = "enable_itr_identification",
        default_value = "false",
        help_heading = "Terminal repeat identification"
    )]
    disable_dtr_trimming: bool,

    /// Minimum length of terminal repeat
    #[clap(
        short = 'l',
        long,
        value_parser,
        default_value = "21",
        help_heading = "Terminal repeat identification"
    )]
    min_length: usize,

    /// Ignore terminal repeats that contain a high proportion of low complexity
    /// sequences
    #[clap(
        long,
        short = 'c',
        value_parser,
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
        value_parser,
        default_value = "false",
        help_heading = "Terminal repeat filtering"
    )]
    ignore_ambiguous: bool,

    /// Maximum fraction of the terminal repeat length that is comprised of
    /// ambiguous bases
    #[clap(
        long,
        value_parser = fraction_in_range,
        default_value = "0.5",
        requires = "ignore_ambiguous",
        help_heading = "Terminal repeat filtering"
    )]
    max_ambiguous_frac: f64,

    /// Retain only the sequences for which terminal repeats were identified
    #[clap(
        short = 'x',
        long,
        value_parser,
        default_value = "false",
        help_heading = "Output"
    )]
    exclude_non_tr_seqs: bool,

    /// Add terminal repeat information to the sequence headers (e.g.,
    /// 'tr=dtr tr_length=55')
    #[clap(
        short = 'a',
        long,
        value_parser,
        default_value = "false",
        help_heading = "Output"
    )]
    include_tr_info: bool,

    /// Disable trimming of terminal repeats from sequences. Can be used with
    /// `--include-tr-info` or `--exclude-non-tr-seqs` to identify and report
    /// sequences with terminal repeats without modifying the sequences
    #[clap(
        short = 't',
        long,
        value_parser,
        default_value = "false",
        help_heading = "Output"
    )]
    disable_trimming: bool,
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
    input: Input,
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
    let reader = match input.is_std() {
        true => parse_fastx_stdin(),
        false => {
            if input.is_empty().unwrap() {
                eprintln!("Error: the input file is empty");
                process::exit(1);
            }
            parse_fastx_file(input.path().to_path_buf())
        }
    };

    let mut reader = match reader {
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    };

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
    let enable_itr_identification = cli.enable_itr_identification;
    let disable_dtr_trimming = cli.disable_dtr_trimming;
    let min_length = cli.min_length;
    let ignore_low_complexity = cli.ignore_low_complexity;
    let max_low_complexity_frac = cli.max_low_complexity_frac;
    let ignore_ambiguous = cli.ignore_ambiguous;
    let max_ambiguous_frac = cli.max_ambiguous_frac;
    let exclude_non_tr_seqs = cli.exclude_non_tr_seqs;
    let include_tr_info = cli.include_tr_info;
    let disable_trimming = cli.disable_trimming;
    for input in cli.input {
        pipeline(
            input,
            enable_itr_identification,
            disable_dtr_trimming,
            min_length,
            ignore_low_complexity,
            max_low_complexity_frac,
            ignore_ambiguous,
            max_ambiguous_frac,
            exclude_non_tr_seqs,
            include_tr_info,
            disable_trimming,
        );
    }
}
