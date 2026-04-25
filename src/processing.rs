use crate::tr::{
    RepeatKind, RepeatOptions, RepeatScratch, RepeatSearchResult, TerminalRepeat,
    find_repeats_with_scratch,
};
use needletail::parser::FastxReader;
use std::io::{self, Write};
use std::str::{Utf8Error, from_utf8};

const FASTA_LINE_WIDTH: usize = 80;

#[derive(Clone, Copy, Debug)]
pub struct FastaOutputOptions {
    pub include_tr_info: bool,
    pub disable_trimming: bool,
}

#[derive(Debug)]
pub enum ProcessError {
    Read(String),
    Write(io::Error),
}

pub fn process_reader(
    mut reader: Box<dyn FastxReader>,
    writer: &mut impl Write,
    repeat_options: RepeatOptions,
    exclude_non_tr_seqs: bool,
    output_options: FastaOutputOptions,
) -> Result<(), ProcessError> {
    let mut repeat_scratch = RepeatScratch::default();

    while let Some(record) = reader.next() {
        let record = record.map_err(|error| ProcessError::Read(error.to_string()))?;
        let sequence = record.seq();

        let repeat =
            find_repeats_with_scratch(sequence.as_ref(), repeat_options, &mut repeat_scratch);

        if exclude_non_tr_seqs && repeat.accepted_repeat().is_none() {
            continue;
        }

        match write_fasta_record(
            writer,
            record.id(),
            sequence.as_ref(),
            repeat,
            output_options,
        ) {
            Ok(()) => {}
            Err(FastaWriteError::Io(error)) => return Err(ProcessError::Write(error)),
            Err(FastaWriteError::Utf8(error)) => {
                eprintln!("Error formatting record: {error}");
            }
        }
    }

    Ok(())
}

fn write_fasta_record(
    writer: &mut impl Write,
    header: &[u8],
    sequence: &[u8],
    repeat: RepeatSearchResult,
    output_options: FastaOutputOptions,
) -> Result<(), FastaWriteError> {
    from_utf8(header)?;
    let accepted_repeat = repeat.accepted_repeat();
    let output_length =
        if let (Some(repeat), false) = (accepted_repeat, output_options.disable_trimming) {
            sequence.len() - repeat.length
        } else {
            sequence.len()
        };

    let sequence = from_utf8(sequence)?;
    let wrapped_sequence = textwrap::fill(&sequence[..output_length], FASTA_LINE_WIDTH);

    writer.write_all(b">")?;
    writer.write_all(header)?;
    if output_options.include_tr_info {
        if let Some(repeat) = accepted_repeat {
            write_repeat_info(writer, repeat)?;
        } else {
            writer.write_all(b" tr=none tr_length=0")?;
        }
    }
    writer.write_all(b"\n")?;

    writer.write_all(wrapped_sequence.as_bytes())?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_repeat_info(writer: &mut impl Write, repeat: TerminalRepeat) -> io::Result<()> {
    write!(
        writer,
        " tr={} tr_length={}",
        repeat_kind_label(repeat.kind),
        repeat.length
    )
}

fn repeat_kind_label(kind: RepeatKind) -> &'static str {
    match kind {
        RepeatKind::Direct => "dtr",
        RepeatKind::Inverted => "itr",
    }
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
