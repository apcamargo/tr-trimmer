use std::collections::VecDeque;

const ENCODING_LOOKUP: [u8; 256] = {
    let mut lookup = [4; 256];

    lookup[b'A' as usize] = 0;
    lookup[b'C' as usize] = 1;
    lookup[b'G' as usize] = 2;
    lookup[b'T' as usize] = 3;
    lookup[b'a' as usize] = 0;
    lookup[b'c' as usize] = 1;
    lookup[b'g' as usize] = 2;
    lookup[b't' as usize] = 3;
    lookup
};

const MASK: u8 = 63;

#[derive(Debug)]
struct PerfectInterval {
    start: usize,
    finish: usize,
    score: usize,
    l: usize,
}

#[derive(Debug)]
struct SymmetricDust<'a> {
    /// `q` in the paper
    sequence: &'a [u8],
    /// The length of the window used by symmetric DUST algorithm
    /// `W` in the paper
    window_size: usize,
    /// 10 times the score threshold used by symmetric DUST algorithm.
    /// `T` in the paper
    score_threshold: usize,
    /// `P` in the paper
    perfect_intervals: VecDeque<PerfectInterval>,
    /// `w` in the paper
    window: VecDeque<usize>,
    // counts in the current window
    cv: [usize; 64],
    cw: [usize; 64],
    // runnings counts
    rv: usize,
    rw: usize,
    /// `L` in the paper
    biggest_num_triplets: usize,
    prefix_limit: usize,
    max_prefix_masked_bases: f64,
    prefix_masked_bases: usize,
    prefix_counted_until: usize,
    prefix_exceeded: bool,
}

impl<'a> SymmetricDust<'a> {
    fn process_prefix_count(
        sequence: &'a [u8],
        window_size: usize,
        score_threshold: usize,
        prefix_len: usize,
        max_prefix_masked_bases: f64,
    ) -> usize {
        let mut obj = SymmetricDust {
            sequence,
            window_size,
            score_threshold,
            perfect_intervals: VecDeque::new(),
            window: VecDeque::new(),
            cv: [0; 64],
            cw: [0; 64],
            rv: 0,
            rw: 0,
            biggest_num_triplets: 0,
            prefix_limit: prefix_len.min(sequence.len()),
            max_prefix_masked_bases,
            prefix_masked_bases: 0,
            prefix_counted_until: 0,
            prefix_exceeded: false,
        };

        obj.inner_process();
        obj.prefix_masked_bases
    }

    fn inner_process(&mut self) {
        // We're going to represent 3 chars in that u8
        let mut triplet: u8 = 0;
        let mut l = 0usize;
        for i in 0..=self.sequence.len() {
            if self.prefix_exceeded {
                break;
            }

            let b = if i < self.sequence.len() {
                ENCODING_LOOKUP[self.sequence[i] as usize]
            } else {
                4
            };

            // A/T/C/G
            if b < 4 {
                l += 1;
                triplet = (triplet << 2 | b) & MASK;

                // We have at least 3 chars, we can look at them
                if l >= 3 {
                    let mut window_start = l.saturating_sub(self.window_size);
                    window_start += i + 1 - l;

                    self.save_masked_regions(window_start);
                    if self.prefix_scan_is_complete(window_start) {
                        break;
                    }
                    self.shift_window(triplet as usize);
                    if self.rw * 10 > self.biggest_num_triplets * self.score_threshold {
                        self.find_perfect(window_start);
                    }
                }
            } else {
                // A `N` resets the sequence
                // When we are there (N or end of seq), we empty the intervals found so far
                let mut window_start = if l > self.window_size - 1 {
                    l - self.window_size + 1
                } else {
                    0
                };
                window_start += i + 1 - l;
                while !self.perfect_intervals.is_empty() {
                    window_start += 1;
                    self.save_masked_regions(window_start);
                    if self.prefix_exceeded || self.prefix_scan_is_complete(window_start) {
                        break;
                    }
                }

                l = 0;
                triplet = 0;
            }
        }
    }

    /// Save all the intervals that are before the `window_start`
    /// This can only insert one result at a time
    fn save_masked_regions(&mut self, window_start: usize) {
        if self.perfect_intervals.is_empty() {
            return;
        }

        let back = self.perfect_intervals.back().unwrap();
        if back.start >= window_start {
            return;
        }
        let (start, finish) = (back.start, back.finish);

        self.record_masked_region(start, finish);

        while let Some(b) = self.perfect_intervals.back() {
            if b.start < window_start {
                self.perfect_intervals.pop_back();
            } else {
                break;
            }
        }
    }

    fn record_masked_region(&mut self, start: usize, finish: usize) {
        if start >= self.prefix_limit {
            return;
        }

        let end = finish.min(self.prefix_limit);
        let count_start = start.max(self.prefix_counted_until);
        if end > count_start {
            self.prefix_masked_bases += end - count_start;
            self.prefix_counted_until = end;
            self.prefix_exceeded = self.prefix_masked_bases as f64 > self.max_prefix_masked_bases;
        }
    }

    fn prefix_scan_is_complete(&self, window_start: usize) -> bool {
        window_start > self.prefix_limit
    }

    /// Add a triplet to the window, shifting all the data to represent the new window
    fn shift_window(&mut self, triplet: usize) {
        let mut s;
        if self.window.len() >= self.window_size - 2 {
            s = self.window.pop_front().unwrap();
            self.cw[s] -= 1;
            self.rw -= self.cw[s];
            if self.biggest_num_triplets > self.window.len() {
                self.biggest_num_triplets -= 1;
                self.cv[s] -= 1;
                self.rv -= self.cv[s];
            }
        }

        self.window.push_back(triplet);
        self.biggest_num_triplets += 1;

        self.rw += self.cw[triplet];
        self.cw[triplet] += 1;
        self.rv += self.cv[triplet];
        self.cv[triplet] += 1;

        if self.cv[triplet] * 10 > 2 * self.score_threshold {
            loop {
                s = self.window[self.window.len() - self.biggest_num_triplets];
                self.biggest_num_triplets -= 1;
                self.cv[s] -= 1;
                self.rv -= self.cv[s];

                if s == triplet {
                    break;
                }
            }
        }
    }

    /// Find all the perfect intervals in the window
    fn find_perfect(&mut self, window_start: usize) {
        let mut c = self.cv;
        let mut r = self.rv;
        let mut max_score = 0;
        let mut max_l = 0;

        for i in (0..self.window.len() - self.biggest_num_triplets).rev() {
            let triplet = self.window[i];
            r += c[triplet];
            c[triplet] += 1;
            let new_score = r;
            let new_l = self.window.len() - i - 1;
            if new_score * 10 > self.score_threshold * new_l {
                let mut insertion_position = 0;
                // Figure out where to insert the new interval
                for (j, interval) in self.perfect_intervals.iter().enumerate() {
                    if interval.start < i + window_start {
                        break;
                    }
                    insertion_position = j + 1;
                    if max_score == 0 || interval.score * max_l > max_score * interval.l {
                        max_score = interval.score;
                        max_l = interval.l;
                    }
                }

                // And insert it
                if max_score == 0 || new_score * max_l >= max_score * new_l {
                    max_score = new_score;
                    max_l = new_l;
                    let new_perf = PerfectInterval {
                        start: i + window_start,
                        // +2 => triplet size (3) - 1
                        finish: self.window.len() + 2 + window_start,
                        score: new_score,
                        l: new_l,
                    };

                    self.perfect_intervals.insert(insertion_position, new_perf);
                }
            }
        }
    }
}

/// Returns how many bases in `sequence[..prefix_len]` are part of low-complexity
/// intervals. Stops once the count exceeds `max_prefix_masked_bases`.
pub fn dustmasker_masked_prefix_bases(
    sequence: &[u8],
    window_size: usize,
    score_threshold: usize,
    prefix_len: usize,
    max_prefix_masked_bases: f64,
) -> usize {
    SymmetricDust::process_prefix_count(
        sequence,
        window_size,
        score_threshold,
        prefix_len,
        max_prefix_masked_bases,
    )
}
