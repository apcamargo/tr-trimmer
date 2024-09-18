use crate::sdust::dustmasker;
use needletail::Sequence;

fn find_dtr(sequence: &[u8], min_length: usize) -> (bool, usize) {
    let seq_len = sequence.len();
    if seq_len < min_length * 2 {
        return (false, 0);
    }
    for length in (min_length..=seq_len / 2).rev() {
        let start = &sequence[..length];
        let end = &sequence[seq_len - length..];
        if start.eq_ignore_ascii_case(end) {
            return (true, length);
        }
    }
    (false, 0)
}

fn find_itr(sequence: &[u8], min_length: usize) -> (bool, usize) {
    let seq_len = sequence.len();
    let rev_complement = sequence.reverse_complement();
    if seq_len < min_length * 2 {
        return (false, 0);
    }
    let start = &sequence[..min_length];
    let end = &rev_complement[..min_length];
    if start.eq(end) {
        let mut i = min_length;
        while i <= seq_len / 2 && sequence[..i].eq_ignore_ascii_case(&rev_complement[..i]) {
            i += 1;
        }
        (true, i - 1)
    } else {
        (false, 0)
    }
}

/// Evaluate the fraction of the TR that is low complexity. Returns false if the
/// fraction of the TR length that is low-complexity exceeds the maximum allowed
/// fraction (`max_lc_frac`).
fn evaluate_tr_complexity(sequence: &[u8], tr_length: usize, max_lc_frac: f64) -> bool {
    let mask = dustmasker(sequence, 32, 30);
    let n_lc_tr: usize = mask
        .iter()
        .take_while(|range| range.start < tr_length)
        .map(|range| range.end.min(tr_length) - range.start)
        .sum();
    (n_lc_tr as f64) / (tr_length as f64) <= max_lc_frac
}

fn evaluate_ambiguous_bases(sequence: &[u8], tr_length: usize, max_ambig_frac: f64) -> bool {
    let norm_sequence = sequence.normalize(false);
    let n_ambig = norm_sequence[..tr_length]
        .iter()
        .filter(|&&base| base == b'N')
        .count();
    (n_ambig as f64) / (tr_length as f64) <= max_ambig_frac
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
    if !disable_dtr_identification {
        let (has_dtr, tr_length) = find_dtr(sequence, min_length);
        if has_dtr || !enable_itr_identification {
            if ignore_low_complexity {
                match evaluate_tr_complexity(sequence, tr_length, max_low_complexity_frac) {
                    true => return (has_dtr, false, tr_length),
                    false => return (false, false, tr_length),
                }
            }
            if ignore_ambiguous {
                match evaluate_ambiguous_bases(sequence, tr_length, max_ambiguous_frac) {
                    true => return (has_dtr, false, tr_length),
                    false => return (false, false, tr_length),
                }
            }
            return (has_dtr, false, tr_length);
        }
    }
    if enable_itr_identification {
        let (has_itr, tr_length) = find_itr(sequence, min_length);
        if ignore_low_complexity {
            match evaluate_tr_complexity(sequence, tr_length, max_low_complexity_frac) {
                true => return (false, has_itr, tr_length),
                false => return (false, false, tr_length),
            }
        }
        if ignore_ambiguous {
            match evaluate_ambiguous_bases(sequence, tr_length, max_ambiguous_frac) {
                true => return (false, has_itr, tr_length),
                false => return (false, false, tr_length),
            }
        }
        return (false, has_itr, tr_length);
    }
    (false, false, 0)
}
