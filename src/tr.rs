use crate::sdust::dustmasker_masked_prefix_bases;
use needletail::Sequence;

const DUST_WINDOW_SIZE: usize = 32;
const DUST_SCORE_THRESHOLD: usize = 30;

#[derive(Clone, Copy, Debug)]
pub struct RepeatOptions {
    pub min_length: usize,
    pub disable_dtr_identification: bool,
    pub enable_itr_identification: bool,
    pub ignore_low_complexity: bool,
    pub max_low_complexity_frac: f64,
    pub ignore_ambiguous: bool,
    pub max_ambiguous_frac: f64,
}

#[derive(Default, Debug)]
pub struct RepeatScratch {
    prefix_function: Vec<usize>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RepeatKind {
    Direct,
    Inverted,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TerminalRepeat {
    pub kind: RepeatKind,
    pub length: usize,
}

impl TerminalRepeat {
    fn new(kind: RepeatKind, length: usize) -> Self {
        Self { kind, length }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RepeatSearchResult {
    Accepted(TerminalRepeat),
    Rejected,
    None,
}

impl RepeatSearchResult {
    pub fn accepted_repeat(self) -> Option<TerminalRepeat> {
        match self {
            Self::Accepted(repeat) => Some(repeat),
            Self::Rejected | Self::None => None,
        }
    }
}

struct TerminalRepeatFinder<'a, 'scratch> {
    sequence: &'a [u8],
    options: RepeatOptions,
    scratch: &'scratch mut RepeatScratch,
}

impl<'a, 'scratch> TerminalRepeatFinder<'a, 'scratch> {
    fn new(
        sequence: &'a [u8],
        options: RepeatOptions,
        scratch: &'scratch mut RepeatScratch,
    ) -> Self {
        Self {
            sequence,
            options,
            scratch,
        }
    }

    fn find_dtr_mut(&mut self) -> (bool, usize) {
        let seq_len = self.sequence.len();
        if self.options.min_length > seq_len / 2 {
            return (false, 0);
        }
        if seq_len == 0 {
            return (true, 0);
        }

        self.scratch.prefix_function.clear();
        self.scratch.prefix_function.resize(seq_len, 0);

        for i in 1..seq_len {
            let mut border = self.scratch.prefix_function[i - 1];
            while border > 0 && !eq_ignore_ascii_case_byte(self.sequence[i], self.sequence[border])
            {
                border = self.scratch.prefix_function[border - 1];
            }
            if eq_ignore_ascii_case_byte(self.sequence[i], self.sequence[border]) {
                border += 1;
            }
            self.scratch.prefix_function[i] = border;
        }

        let mut length = self.scratch.prefix_function[seq_len - 1];
        while length > seq_len / 2 {
            length = self.scratch.prefix_function[length - 1];
        }

        if length >= self.options.min_length {
            (true, length)
        } else {
            (false, 0)
        }
    }

    fn find_itr(&self) -> (bool, usize) {
        let seq_len = self.sequence.len();
        if self.options.min_length > seq_len / 2 {
            return (false, 0);
        }

        let max_length = seq_len / 2;
        let mut length = 0;
        while length < max_length
            && is_itr_pair(self.sequence[length], self.sequence[seq_len - length - 1])
        {
            length += 1;
        }

        if length >= self.options.min_length {
            (true, length)
        } else {
            (false, 0)
        }
    }

    /// Evaluate the fraction of the TR that is low complexity. Returns false if
    /// the fraction of the TR length that is low-complexity exceeds the maximum
    /// allowed fraction (`max_lc_frac`).
    fn is_valid_complexity(&self, tr_length: usize) -> bool {
        if tr_length == 0 {
            return false;
        }

        // If the sequence is longer than 50 * tr_length, dustmasker will
        // process the first 50 * tr_length bases. Otherwise, it will process
        // the entire sequence.
        let max_scan_len = 50 * tr_length;
        let sequence = if self.sequence.len() > max_scan_len {
            &self.sequence[..max_scan_len]
        } else {
            self.sequence
        };
        let max_low_complexity_bases = self.options.max_low_complexity_frac * tr_length as f64;
        let n_lc_tr = dustmasker_masked_prefix_bases(
            sequence,
            DUST_WINDOW_SIZE,
            DUST_SCORE_THRESHOLD,
            tr_length,
            max_low_complexity_bases,
        );
        (n_lc_tr as f64) / (tr_length as f64) <= self.options.max_low_complexity_frac
    }

    fn is_valid_ambiguous_bases(&self, tr_length: usize) -> bool {
        if tr_length == 0 {
            return false;
        }

        if self
            .sequence
            .iter()
            .take(tr_length)
            .any(is_ascii_whitespace)
        {
            let norm_sequence = self.sequence.normalize(false);
            let n_ambig = norm_sequence[..tr_length]
                .iter()
                .filter(|&&base| base == b'N')
                .count();
            return (n_ambig as f64) / (tr_length as f64) <= self.options.max_ambiguous_frac;
        }

        let max_ambiguous_bases = self.options.max_ambiguous_frac * tr_length as f64;
        let mut n_ambig = 0;
        for &base in self.sequence.iter().take(tr_length) {
            if normalizes_to_ambiguous_base(base) {
                n_ambig += 1;
                if n_ambig as f64 > max_ambiguous_bases {
                    return false;
                }
            }
        }
        true
    }

    fn validate_repeat(
        &self,
        tr_length: usize,
        check_complexity: bool,
        check_ambiguous: bool,
    ) -> bool {
        if check_complexity && !self.is_valid_complexity(tr_length) {
            return false;
        }
        if check_ambiguous && !self.is_valid_ambiguous_bases(tr_length) {
            return false;
        }
        true
    }
}

fn eq_ignore_ascii_case_byte(left: u8, right: u8) -> bool {
    left.eq_ignore_ascii_case(&right)
}

fn is_itr_pair(prefix_base: u8, suffix_base: u8) -> bool {
    eq_ignore_ascii_case_byte(prefix_base, complement_base(suffix_base))
}

fn complement_base(base: u8) -> u8 {
    match base {
        b'a' => b't',
        b'A' => b'T',
        b'c' => b'g',
        b'C' => b'G',
        b'g' => b'c',
        b'G' => b'C',
        b't' => b'a',
        b'T' => b'A',
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        b's' => b's',
        b'w' => b'w',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'S' => b'S',
        b'W' => b'W',
        other => other,
    }
}

fn normalizes_to_ambiguous_base(base: u8) -> bool {
    !matches!(
        base,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'a'
            | b'c'
            | b'g'
            | b't'
            | b'u'
            | b'U'
            | b'-'
            | b'.'
            | b'~'
            | b' '
            | b'\t'
            | b'\r'
            | b'\n'
    )
}

fn is_ascii_whitespace(base: &u8) -> bool {
    matches!(*base, b' ' | b'\t' | b'\r' | b'\n' | 0x0b | 0x0c)
}

fn repeat_search_result(kind: RepeatKind, length: usize, is_valid: bool) -> RepeatSearchResult {
    if is_valid {
        RepeatSearchResult::Accepted(TerminalRepeat::new(kind, length))
    } else {
        RepeatSearchResult::Rejected
    }
}

pub fn find_repeats_with_scratch(
    sequence: &[u8],
    options: RepeatOptions,
    scratch: &mut RepeatScratch,
) -> RepeatSearchResult {
    let mut finder = TerminalRepeatFinder::new(sequence, options, scratch);
    if !options.disable_dtr_identification {
        let (has_dtr, tr_length) = finder.find_dtr_mut();
        if has_dtr || !options.enable_itr_identification {
            let is_valid = has_dtr
                && finder.validate_repeat(
                    tr_length,
                    options.ignore_low_complexity,
                    options.ignore_ambiguous,
                );
            if has_dtr {
                return repeat_search_result(RepeatKind::Direct, tr_length, is_valid);
            }
            return RepeatSearchResult::None;
        }
    }
    if options.enable_itr_identification {
        let (has_itr, tr_length) = finder.find_itr();
        if has_itr {
            let is_valid = finder.validate_repeat(
                tr_length,
                options.ignore_low_complexity,
                options.ignore_ambiguous,
            );
            return repeat_search_result(RepeatKind::Inverted, tr_length, is_valid);
        }
    }
    RepeatSearchResult::None
}
