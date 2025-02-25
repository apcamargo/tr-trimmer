use crate::sdust::dustmasker;
use needletail::Sequence;

struct TerminalRepeatFinder<'a> {
    sequence: &'a [u8],
    min_length: usize,
    max_lc_frac: f64,
    max_ambig_frac: f64,
}

impl<'a> TerminalRepeatFinder<'a> {
    fn new(sequence: &'a [u8], min_length: usize, max_lc_frac: f64, max_ambig_frac: f64) -> Self {
        Self {
            sequence,
            min_length,
            max_lc_frac,
            max_ambig_frac,
        }
    }

    fn find_dtr(&self) -> (bool, usize) {
        let seq_len = self.sequence.len();
        if seq_len < self.min_length * 2 {
            return (false, 0);
        }
        for length in (self.min_length..=seq_len / 2).rev() {
            let start = &self.sequence[..length];
            let end = &self.sequence[seq_len - length..];
            if start.eq_ignore_ascii_case(end) {
                return (true, length);
            }
        }
        (false, 0)
    }

    fn find_itr(&self) -> (bool, usize) {
        let seq_len = self.sequence.len();
        let rev_complement = self.sequence.reverse_complement();
        if seq_len < self.min_length * 2 {
            return (false, 0);
        }
        let start = &self.sequence[..self.min_length];
        let end = &rev_complement[..self.min_length];
        if !start.eq_ignore_ascii_case(end) {
            return (false, 0);
        }
        let mut length = self.min_length;
        while length <= seq_len / 2
            && self.sequence[..length].eq_ignore_ascii_case(&rev_complement[..length])
        {
            length += 1;
        }
        (true, length - 1)
    }

    /// Evaluate the fraction of the TR that is low complexity. Returns false if
    /// the fraction of the TR length that is low-complexity exceeds the maximum
    /// allowed fraction (`max_lc_frac`).
    fn is_valid_complexity(&self, tr_length: usize) -> bool {
        // If the sequence is longer than 50 * tr_length, dustmasker will
        // process the first 50 * tr_length bases. Otherwise, it will process
        // the entire sequence.
        let mask = if self.sequence.len() > 50 * tr_length {
            dustmasker(&self.sequence[..50 * tr_length], 32, 30)
        } else {
            dustmasker(&self.sequence, 32, 30)
        };
        let n_lc_tr: usize = mask
            .iter()
            .take_while(|range| range.start < tr_length)
            .map(|range| range.end.min(tr_length) - range.start)
            .sum();
        (n_lc_tr as f64) / (tr_length as f64) <= self.max_lc_frac
    }

    fn is_valid_ambiguous_bases(&self, tr_length: usize) -> bool {
        let norm_sequence = self.sequence.normalize(false);
        let n_ambig = norm_sequence[..tr_length]
            .iter()
            .filter(|&&base| base == b'N')
            .count();
        (n_ambig as f64) / (tr_length as f64) <= self.max_ambig_frac
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

pub fn find_repeats(
    sequence: &[u8],
    min_length: usize,
    disable_dtr_identification: bool,
    enable_itr_identification: bool,
    ignore_low_complexity: bool,
    max_low_complexity_frac: f64,
    ignore_ambiguous: bool,
    max_ambiguous_frac: f64,
) -> (bool, bool, usize) {
    let finder = TerminalRepeatFinder::new(
        sequence,
        min_length,
        max_low_complexity_frac,
        max_ambiguous_frac,
    );
    if !disable_dtr_identification {
        let (has_dtr, tr_length) = finder.find_dtr();
        if has_dtr || !enable_itr_identification {
            let is_valid =
                finder.validate_repeat(tr_length, ignore_low_complexity, ignore_ambiguous);
            return (is_valid && has_dtr, false, tr_length);
        }
    }
    if enable_itr_identification {
        let (has_itr, tr_length) = finder.find_itr();
        if has_itr {
            let is_valid =
                finder.validate_repeat(tr_length, ignore_low_complexity, ignore_ambiguous);
            return (false, is_valid && has_itr, tr_length);
        }
    }
    (false, false, 0)
}
