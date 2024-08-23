use std::collections::VecDeque;

pub(crate) struct LogThinStats {
    n: usize,
    x_sum: f64,
    var_sum: f64,
    mean: f64,
    variance: f64,
    auto_corrs: AutoCorrs,
}

const SHORT_LAG_BUF_LENGTH: usize = 16;
const N_BIN_LEVELS: usize = 10;

struct ShortLags {
    corr_sums: [f64; SHORT_LAG_BUF_LENGTH],
    x_diffs: VecDeque<f64>,
}

struct Bins {
   bins: [[f64; 2]; N_BIN_LEVELS],
}
struct LongLags {
    bins_for_odd: Bins,
    bins_for_even: Bins,
}
struct AutoCorrs {
    short_lags: ShortLags,
    long_lags: LongLags
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
        self.n += 1;
        self.x_sum += x;
        self.var_sum += x.powi(2);
        let x_mean_previous = self.mean;
        self.mean = self.x_sum / (self.n as f64);
        let x_diff = x - self.mean;
        self.var_sum += (x - x_mean_previous) * x_diff;
        self.variance = self.var_sum / (self.n as f64);
    }
    pub(crate) fn mean(&self) -> f64 { self.mean }
    pub(crate) fn variance(&self) -> f64 { self.variance }
}

impl ShortLags {
    pub fn new() -> ShortLags {
        let corr_sums = [0.0; SHORT_LAG_BUF_LENGTH];
        let x_diffs = VecDeque::with_capacity(SHORT_LAG_BUF_LENGTH);
        ShortLags { corr_sums, x_diffs }
    }
    pub fn add_x_diff(&mut self, x_diff: f64) {
        self.x_diffs.push_back(x_diff);
    }
}

impl Bins {
    pub fn new() -> Bins {
        let bins = [[0.0; 2]; N_BIN_LEVELS];
        Bins { bins }
    }
}

impl LongLags {
    pub fn new() -> LongLags {
        let bins_for_odd = Bins::new();
        let bins_for_even = Bins::new();
        LongLags { bins_for_odd, bins_for_even }
    }
    pub fn add_x_diff(&mut self, x_diff: f64) {
    }
}

impl AutoCorrs {
    pub(crate) fn new() -> AutoCorrs {
        let short_lags = ShortLags::new();
        let long_lags = LongLags::new();
        AutoCorrs { short_lags, long_lags }
    }
    pub(crate) fn add_x_diff(&mut self, x_diff: f64) {
        self.short_lags.add_x_diff(x_diff);
    }
}
