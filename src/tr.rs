use crate::sdust::dustmasker_masked_prefix_bases;
use needletail::Sequence;

const DUST_WINDOW_SIZE: usize = 32;
const DUST_SCORE_THRESHOLD: usize = 30;

#[derive(Clone, Copy, Debug)]
pub(crate) struct RepeatOptions {
    pub(crate) min_length: usize,
    pub(crate) disable_dtr_identification: bool,
    pub(crate) enable_itr_identification: bool,
    pub(crate) ignore_low_complexity: bool,
    pub(crate) max_low_complexity_frac: f64,
    pub(crate) ignore_ambiguous: bool,
    pub(crate) max_ambiguous_frac: f64,
}

#[derive(Default, Debug)]
pub(crate) struct RepeatScratch {
    prefix_function: Vec<usize>,
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

pub(crate) fn find_repeats_with_scratch(
    sequence: &[u8],
    options: RepeatOptions,
    scratch: &mut RepeatScratch,
) -> (bool, bool, usize) {
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
            return (is_valid, false, tr_length);
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
            return (false, is_valid, tr_length);
        }
    }
    (false, false, 0)
}

#[cfg(test)]
fn find_repeats(sequence: &[u8], options: RepeatOptions) -> (bool, bool, usize) {
    let mut scratch = RepeatScratch::default();
    find_repeats_with_scratch(sequence, options, &mut scratch)
}

#[cfg(test)]
mod tests {
    use super::{
        RepeatOptions, RepeatScratch, TerminalRepeatFinder, find_repeats, find_repeats_with_scratch,
    };
    use needletail::Sequence;

    fn options(min_length: usize) -> RepeatOptions {
        RepeatOptions {
            min_length,
            disable_dtr_identification: false,
            enable_itr_identification: false,
            ignore_low_complexity: false,
            max_low_complexity_frac: 0.5,
            ignore_ambiguous: false,
            max_ambiguous_frac: 0.0,
        }
    }

    fn optimized_dtr(sequence: &[u8], min_length: usize) -> (bool, usize) {
        let mut scratch = RepeatScratch::default();
        let mut finder = TerminalRepeatFinder::new(sequence, options(min_length), &mut scratch);
        finder.find_dtr_mut()
    }

    fn optimized_itr(sequence: &[u8], min_length: usize) -> (bool, usize) {
        let mut scratch = RepeatScratch::default();
        let finder = TerminalRepeatFinder::new(sequence, options(min_length), &mut scratch);
        finder.find_itr()
    }

    fn reference_dtr(sequence: &[u8], min_length: usize) -> (bool, usize) {
        if sequence.len() < min_length.saturating_mul(2) {
            return (false, 0);
        }

        for length in (min_length..=sequence.len() / 2).rev() {
            if sequence[..length].eq_ignore_ascii_case(&sequence[sequence.len() - length..]) {
                return (true, length);
            }
        }
        (false, 0)
    }

    fn reference_itr(sequence: &[u8], min_length: usize) -> (bool, usize) {
        if sequence.len() < min_length.saturating_mul(2) {
            return (false, 0);
        }

        let rev_complement = sequence.reverse_complement();
        if !sequence[..min_length].eq_ignore_ascii_case(&rev_complement[..min_length]) {
            return (false, 0);
        }

        let mut length = min_length;
        while length <= sequence.len() / 2
            && sequence[..length].eq_ignore_ascii_case(&rev_complement[..length])
        {
            length += 1;
        }
        (true, length - 1)
    }

    fn repeated_prefix(pattern: &[u8], length: usize) -> Vec<u8> {
        pattern.iter().copied().cycle().take(length).collect()
    }

    fn deterministic_sequence(seed: u64, length: usize, alphabet: &[u8]) -> Vec<u8> {
        let mut state = seed.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut sequence = Vec::with_capacity(length);
        for _ in 0..length {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            sequence.push(alphabet[(state as usize >> 16) % alphabet.len()]);
        }
        sequence
    }

    fn reference_sequences() -> Vec<Vec<u8>> {
        let mut sequences = vec![
            Vec::new(),
            b"A".to_vec(),
            b"AAAAA".to_vec(),
            b"ACGTACGT".to_vec(),
            b"acgtNNACGT".to_vec(),
            b"ACGNttACGN".to_vec(),
            b"RYKMACCKMRY".to_vec(),
            b"ACGUACGUxxACGTACGT".to_vec(),
            b"ACG-ACG-xxACG-ACG-".to_vec(),
            b"ACG.ACG~xxACG.ACG~".to_vec(),
            b"ACGT ACGTxxACGT ACGT".to_vec(),
            b"ACGT\tACGTxxACGT\tACGT".to_vec(),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
            b"ATATATATATATCGATATATATATAT".to_vec(),
        ];

        for length in [
            0, 1, 2, 3, 4, 5, 8, 20, 21, 22, 31, 32, 33, 63, 64, 65, 79, 80, 81, 127, 128, 129,
            255, 256,
        ] {
            let repeat = repeated_prefix(b"ACGTRYKMBDHVSWN", length);
            sequences.push([repeat.as_slice(), b"NNGATTACA", repeat.as_slice()].concat());
            sequences.push(
                [
                    repeat.to_ascii_lowercase().as_slice(),
                    b"NNGATTACA",
                    repeat.to_ascii_uppercase().as_slice(),
                ]
                .concat(),
            );
            sequences.push(
                [
                    repeat.as_slice(),
                    b"CCGGTA",
                    repeat.reverse_complement().as_slice(),
                ]
                .concat(),
            );
            sequences.push(
                [
                    repeated_prefix(b"AC-GT", length).as_slice(),
                    b"xx",
                    repeated_prefix(b"AC-GT", length).as_slice(),
                ]
                .concat(),
            );
        }

        for seed in 0..64 {
            sequences.push(deterministic_sequence(
                seed,
                1 + seed as usize * 7,
                b"ACGTN",
            ));
            sequences.push(deterministic_sequence(
                seed + 100,
                2 + seed as usize * 5,
                b"ACGTRYKMBDHVSWNacgtnu.-~",
            ));
        }

        sequences
    }

    #[test]
    fn optimized_dtr_matches_naive_reference_for_complex_sequences() {
        for sequence in reference_sequences() {
            for min_length in [
                0, 1, 2, 3, 4, 5, 20, 21, 22, 31, 32, 33, 63, 64, 65, 79, 80, 81, 127, 128, 129,
                255, 256,
            ] {
                assert_eq!(
                    optimized_dtr(&sequence, min_length),
                    reference_dtr(&sequence, min_length),
                    "sequence={} min_length={min_length}",
                    String::from_utf8_lossy(&sequence).escape_debug()
                );
            }
        }
    }

    #[test]
    fn optimized_itr_matches_reverse_complement_reference_for_complex_sequences() {
        for sequence in reference_sequences() {
            for min_length in [
                0, 1, 2, 3, 4, 5, 20, 21, 22, 31, 32, 33, 63, 64, 65, 79, 80, 81, 127, 128, 129,
                255, 256,
            ] {
                assert_eq!(
                    optimized_itr(&sequence, min_length),
                    reference_itr(&sequence, min_length),
                    "sequence={} min_length={min_length}",
                    String::from_utf8_lossy(&sequence).escape_debug()
                );
            }
        }
    }

    #[test]
    fn finds_longest_case_insensitive_dtr() {
        assert_eq!(find_repeats(b"acgtNNACGT", options(4)), (true, false, 4));
    }

    #[test]
    fn rejects_dtr_longer_than_half_sequence() {
        assert_eq!(find_repeats(b"AAAAA", options(3)), (false, false, 0));
    }

    #[test]
    fn finds_itr_without_allocating_reverse_complement() {
        let mut options = options(4);
        options.disable_dtr_identification = true;
        options.enable_itr_identification = true;
        let mut scratch = RepeatScratch::default();
        assert_eq!(
            find_repeats_with_scratch(b"ACGTCCACGT", options, &mut scratch),
            (false, true, 4)
        );
    }

    #[test]
    fn preserves_iupac_complement_behavior_for_itr() {
        let mut options = options(4);
        options.disable_dtr_identification = true;
        options.enable_itr_identification = true;
        let mut scratch = RepeatScratch::default();
        assert_eq!(
            find_repeats_with_scratch(b"RYKMACCKMRY", options, &mut scratch),
            (false, true, 4)
        );
    }

    #[test]
    fn invalid_dtr_filter_does_not_fall_back_to_itr() {
        let mut options = options(4);
        options.enable_itr_identification = true;
        options.ignore_ambiguous = true;
        options.max_ambiguous_frac = 0.0;
        let mut scratch = RepeatScratch::default();
        assert_eq!(
            find_repeats_with_scratch(b"ACGNttACGN", options, &mut scratch),
            (false, false, 4)
        );
    }

    #[test]
    fn zero_length_repeat_is_invalid_when_filters_are_enabled() {
        let mut low_complexity_options = options(0);
        low_complexity_options.ignore_low_complexity = true;
        let mut ambiguous_options = options(0);
        ambiguous_options.ignore_ambiguous = true;

        assert_eq!(find_repeats(b"ACCGGTAT", options(0)), (true, false, 0));
        assert_eq!(
            find_repeats(b"ACCGGTAT", low_complexity_options),
            (false, false, 0)
        );
        assert_eq!(
            find_repeats(b"ACCGGTAT", ambiguous_options),
            (false, false, 0)
        );
    }
}
