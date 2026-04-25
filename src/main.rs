mod processing;
mod sdust;
mod tr;

use crate::processing::{FastaOutputOptions, ProcessError, process_reader};
use crate::tr::RepeatOptions;
use clap::{
    CommandFactory, Parser,
    builder::styling::{AnsiColor, Style, Styles},
};

use clio::Input;
use needletail::{parse_fastx_reader, parser::FastxReader};
use std::io::{self, BufWriter, IsTerminal, Write};
use std::ops::RangeInclusive;
use std::process;

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
    /// (requires --enable-itr-identification)
    #[clap(
        short = 'd',
        long = "disable-dtr-trimming",
        requires = "enable_itr_identification",
        default_value = "false",
        help_heading = "Terminal repeat identification"
    )]
    disable_dtr_identification: bool,

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

fn print_help_or_exit(code: i32) -> ! {
    if let Err(error) = Cli::command().print_help() {
        handle_output_error(error);
    }
    process::exit(code);
}

fn print_help_after_flush_or_exit(writer: &mut impl Write, code: i32) -> ! {
    if let Err(error) = Cli::command().print_help() {
        handle_output_error(error);
    }
    exit_after_flush(writer, code);
}

fn handle_process_error(error: ProcessError, writer: &mut impl Write) -> ! {
    match error {
        ProcessError::Read(error) => {
            eprintln!("Error: {error}");
            exit_after_flush(writer, 1);
        }
        ProcessError::Write(error) => handle_output_error(error),
    }
}

fn main() {
    let cli = Cli::parse();
    let input_count = cli.input.len();

    // If it's an interactive session with no data piped to stdin and files provided,
    // show help and exit
    if input_count == 1 && cli.input[0].is_std() && io::stdin().is_terminal() {
        print_help_or_exit(0);
    }

    let stdout = io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    let repeat_options = RepeatOptions {
        min_length: cli.min_length,
        disable_dtr_identification: cli.disable_dtr_identification,
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
                        print_help_after_flush_or_exit(&mut writer, 0);
                    }
                    // If stdin is invalid but there are other inputs, skip it
                    continue;
                }
                // If the error is from a file input, report and exit
                eprintln!("Error: failed to create reader for {input_display}: {error_msg}");
                exit_after_flush(&mut writer, 1);
            }
        };

        if let Err(error) = process_reader(
            reader,
            &mut writer,
            repeat_options,
            cli.exclude_non_tr_seqs,
            output_options,
        ) {
            handle_process_error(error, &mut writer);
        }
    }

    flush_writer_or_exit(&mut writer);
}
