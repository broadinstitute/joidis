use std::collections::VecDeque;

pub(crate) struct LogThinStats {
    n: usize,
    x_sum: f64,
    var_sum: f64,
    mean: f64,
    variance: f64,
    auto_corrs: AutoCorrs,
}

const SHORT_LAG_BUF_LENGTH: usize = 7;
const N_BIN2_LEVELS: usize = 10;
const N_BIN3_LEVELS: usize = 7;

struct ShortLags {
    corr_sums: [f64; SHORT_LAG_BUF_LENGTH],
    x_diffs: VecDeque<f64>,
}

struct Bins2 {
    bins: [[f64; 2]; N_BIN2_LEVELS],
    corr_sums: [f64; N_BIN2_LEVELS],
}

struct Bins3 {
    bins: [[f64; 3]; N_BIN3_LEVELS],
    corr_sums: [f64; N_BIN3_LEVELS],
}

struct LongLags {
    bins_for_odd: Bins3,
    bins_for_even: Bins2,
}
struct AutoCorrs {
    short_lags: ShortLags,
    long_lags: LongLags,
}

impl LogThinStats {
    pub(crate) fn new() -> LogThinStats {
        let n: usize = 0;
        let x_sum: f64 = 0.0;
        let var_sum: f64 = 0.0;
        let mean: f64 = 0.0;
        let variance: f64 = 0.0;
        let auto_corrs = AutoCorrs::new();
        LogThinStats { n, x_sum, var_sum, mean, variance, auto_corrs }
    }
    pub(crate) fn add(&mut self, x: f64) {
        let n_old = self.n;
        self.n += 1;
        self.x_sum += x;
        self.var_sum += x.powi(2);
        let x_mean_previous = self.mean;
        self.mean = self.x_sum / (self.n as f64);
        let x_diff = x - self.mean;
        self.var_sum += (x - x_mean_previous) * x_diff;
        self.variance = self.var_sum / (self.n as f64);
        self.auto_corrs.add_x_diff(n_old, x_diff);
    }
    pub(crate) fn mean(&self) -> f64 { self.mean }
    pub(crate) fn variance(&self) -> f64 { self.variance }
    pub(crate) fn auto_corr(&self, lag: usize) -> f64 { self.auto_corrs.auto_corr(self.n, lag) }
}

impl ShortLags {
    pub fn new() -> ShortLags {
        let corr_sums = [0.0; SHORT_LAG_BUF_LENGTH];
        let x_diffs = VecDeque::with_capacity(SHORT_LAG_BUF_LENGTH);
        ShortLags { corr_sums, x_diffs }
    }
    pub fn add_x_diff(&mut self, x_diff: f64) {
        for (j, x_diff_j) in self.x_diffs.iter().enumerate() {
            self.corr_sums[j] += x_diff * x_diff_j;
        }
        if self.x_diffs.len() == SHORT_LAG_BUF_LENGTH {
            self.x_diffs.pop_back();
        }
        self.x_diffs.push_front(x_diff);
    }
}

impl Bins2 {
    fn new() -> Bins2 {
        let bins = [[0.0; 2]; N_BIN2_LEVELS];
        let corr_sums = [0.0; N_BIN2_LEVELS];
        Bins2 { bins, corr_sums }
    }
    fn add_x_diff(&mut self, i: usize, x_diff: f64) {
        let mut i = i;
        let mut sum = x_diff;
        for bin_level in 0..N_BIN2_LEVELS {
            let mut level_bins = self.bins[bin_level];
            let i_bin = i % 2;
            if i % 2 == 0 {
                let sum_new = level_bins[0] + level_bins[1];
                level_bins[i_bin] = sum;
                sum = sum_new;
                i = i / 2;
            } else {
                level_bins[i_bin] = sum;
                break;
            }
        }
    }
}

impl Bins3 {
    fn new() -> Bins3 {
        let bins = [[0.0; 3]; N_BIN3_LEVELS];
        let corr_sums = [0.0; N_BIN3_LEVELS];
        Bins3 { bins, corr_sums }
    }
    fn add_x_diff(&mut self, i: usize, x_diff: f64) {
        let mut i = i;
        let mut sum = x_diff;
        for bin_level in 0..N_BIN3_LEVELS {
            let mut level_bins = self.bins[bin_level];
            let i_bin = i % 3;
            if i % 3 == 0 {
                let sum_new = level_bins[0] + level_bins[1] + level_bins[2];
                level_bins[i_bin] = sum;
                sum = sum_new;
                i = i / 3;
            } else {
                level_bins[i_bin] = sum;
                break;
            }
        }
    }
}

impl LongLags {
    fn new() -> LongLags {
        let bins_for_odd = Bins3::new();
        let bins_for_even = Bins2::new();
        LongLags { bins_for_odd, bins_for_even }
    }
    pub fn add_x_diff(&mut self, x_diff: f64) {}
}

impl AutoCorrs {
    fn new() -> AutoCorrs {
        let short_lags = ShortLags::new();
        let long_lags = LongLags::new();
        AutoCorrs { short_lags, long_lags }
    }
    fn add_x_diff(&mut self, i: usize, x_diff: f64) {
        self.short_lags.add_x_diff(x_diff);
    }
    fn auto_corr(&self, n: usize, lag: usize) -> f64 {
        if lag <= SHORT_LAG_BUF_LENGTH {
            self.short_lags.corr_sums[lag] / ((n + lag) as f64)
        } else {
            todo!()
        }
    }
}
