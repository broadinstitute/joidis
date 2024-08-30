use std::collections::VecDeque;

pub(crate) struct PolationStats {
    n: usize,
    x_sum: f64,
    var_sum: f64,
    mean: f64,
    variance: f64,
    covariances: Covariances,
}

const SHORT_LAG_BUF_LENGTH: usize = 32;
const N_BIN2_LEVELS: usize = 10;
const N_BIN3_LEVELS: usize = 7;

struct ShortLags {
    cov_sums: [f64; SHORT_LAG_BUF_LENGTH],
    x_diffs: VecDeque<f64>,
}

struct Bins2 {
    bins: [[f64; 2]; N_BIN2_LEVELS],
    cov_sums: [f64; N_BIN2_LEVELS],
}

struct Bins3 {
    bins: [[f64; 3]; N_BIN3_LEVELS],
    cov_sums: [f64; N_BIN3_LEVELS],
}

struct LongLags {
    bins_for_odd: Bins3,
    bins_for_even: Bins2,
}
struct Covariances {
    short_lags: ShortLags,
    long_lags: LongLags,
}

enum Anchors {
    Single(usize),
    Double(usize, usize),
}

impl PolationStats {
    pub(crate) fn new() -> PolationStats {
        let n: usize = 0;
        let x_sum: f64 = 0.0;
        let var_sum: f64 = 0.0;
        let mean: f64 = 0.0;
        let variance: f64 = 0.0;
        let auto_corrs = Covariances::new();
        PolationStats { n, x_sum, var_sum, mean, variance, covariances: auto_corrs }
    }
    pub(crate) fn add(&mut self, x: f64) {
        let n_old = self.n;
        self.n += 1;
        self.x_sum += x;
        let x_mean_previous = self.mean;
        self.mean = self.x_sum / (self.n as f64);
        let x_diff = x - self.mean;
        self.var_sum += (x - x_mean_previous) * x_diff;
        self.variance = self.var_sum / ((self.n - 1) as f64);
        self.covariances.add_x_diff(n_old, x_diff);
    }
    pub(crate) fn mean(&self) -> f64 { self.mean }
    pub(crate) fn variance(&self) -> f64 { self.variance }
    fn auto_corr_short_lags(&self, lag: usize) -> f64 {
        if lag == 0 {
            1.0
        } else {
            self.covariances.covariance(self.n, lag) / self.variance
        }
    }
    pub(crate) fn auto_corr(&self, lag: usize) -> f64 {
        if lag <= SHORT_LAG_BUF_LENGTH {
            self.auto_corr_short_lags(lag)
        } else {
            let corr_max = self.auto_corr_short_lags(SHORT_LAG_BUF_LENGTH);
            let corr_frac = self.auto_corr_short_lags(lag % SHORT_LAG_BUF_LENGTH);
            corr_max.powi((lag / SHORT_LAG_BUF_LENGTH) as i32) * corr_frac
        }
    }
    pub(crate) fn auto_corr_sum(&self) -> f64 {
        let mut sum = 0.0;
        for lag in 0..SHORT_LAG_BUF_LENGTH {
            sum += self.auto_corr(lag);
        }
        let auto_corr_max = self.auto_corr(SHORT_LAG_BUF_LENGTH);
        sum *= 1.0 / (1.0 - auto_corr_max);
        sum
    }
    pub(crate) fn sum<'a, I>(stats_iter: &mut I) -> PolationStats
    where I: Iterator<Item=&'a PolationStats> {
        let mut n: usize = 0;
        let mut x_sum: f64 = 0.0;
        let mut var_sum: f64 = 0.0;
        let mut sum = PolationStats::new();
        let mut covariances = Covariances::new();
        for stat in stats_iter {
            n += stat.n;
            x_sum += stat.x_sum;
            var_sum += stat.var_sum;
            covariances.add_covariance(&stat.covariances);
        }
        let mean = x_sum / (n as f64);
        let variance = var_sum / ((n - 1) as f64);
        sum.mean = sum.x_sum / (sum.n as f64);
        sum.variance = sum.var_sum / ((sum.n - 1) as f64);
        PolationStats { n, x_sum, var_sum, mean, variance, covariances }
    }
}

impl ShortLags {
    pub fn new() -> ShortLags {
        let cov_sums = [0.0; SHORT_LAG_BUF_LENGTH];
        let x_diffs = VecDeque::with_capacity(SHORT_LAG_BUF_LENGTH);
        ShortLags { cov_sums, x_diffs }
    }
    pub fn add_x_diff(&mut self, x_diff: f64) {
        for (j, x_diff_j) in self.x_diffs.iter().enumerate() {
            self.cov_sums[j] += x_diff * x_diff_j;
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
        let cov_sums = [0.0; N_BIN2_LEVELS];
        Bins2 { bins, cov_sums }
    }
    fn add_x_diff(&mut self, i: usize, x_diff: f64) {
        let mut i = i;
        let mut sum = x_diff;
        for bin_level in 0..N_BIN2_LEVELS {
            let mut level_bins = self.bins[bin_level];
            let i_bin = i % 2;
            if i % 2 == 0 {
                let sum_new = level_bins[0] + level_bins[1];
                let sum_old = level_bins[i_bin];
                level_bins[i_bin] = sum;
                self.cov_sums[bin_level] += x_diff * (sum - sum_old);
                sum = sum_new;
                i /= 2;
            } else {
                level_bins[i_bin] = sum;
                break;
            }
        }
    }
    fn anchors_for(lag: usize) -> Anchors {
        for level1 in 0..(N_BIN2_LEVELS - 1) {
            let level2 = level1 + 1;
            let lag2 = 2_usize.pow(level2 as u32);
            if lag == lag2 {
                return Anchors::Single(level2);
            } else if lag < lag2 {
                return Anchors::Double(level1, level2);
            }
        }
        Anchors::Double(2_usize.pow((N_BIN2_LEVELS - 1) as u32),
                        2_usize.pow(N_BIN2_LEVELS as u32))
    }
}

impl Bins3 {
    fn new() -> Bins3 {
        let bins = [[0.0; 3]; N_BIN3_LEVELS];
        let cov_sums = [0.0; N_BIN3_LEVELS];
        Bins3 { bins, cov_sums }
    }
    fn add_x_diff(&mut self, i: usize, x_diff: f64) {
        let mut i = i;
        let mut sum = x_diff;
        for bin_level in 0..N_BIN3_LEVELS {
            let mut level_bins = self.bins[bin_level];
            let i_bin = i % 3;
            if i % 3 == 0 {
                let sum_new = level_bins[0] + level_bins[1] + level_bins[2];
                let sum_old = level_bins[i_bin];
                level_bins[i_bin] = sum;
                self.cov_sums[bin_level] += x_diff * (sum - sum_old);
                sum = sum_new;
                i /= 3;
            } else {
                level_bins[i_bin] = sum;
                break;
            }
        }
    }
    fn anchors_for(lag: usize) -> Anchors {
        for level1 in 0..(N_BIN3_LEVELS - 1) {
            let level2 = level1 + 1;
            let lag2 = 3_usize.pow(level2 as u32);
            if lag == lag2 {
                return Anchors::Single(level2);
            } else if lag < lag2 {
                return Anchors::Double(level1, level2);
            }
        }
        Anchors::Double(3_usize.pow((N_BIN3_LEVELS - 1) as u32),
                        3_usize.pow(N_BIN3_LEVELS as u32))
    }
}

impl LongLags {
    fn new() -> LongLags {
        let bins_for_odd = Bins3::new();
        let bins_for_even = Bins2::new();
        LongLags { bins_for_odd, bins_for_even }
    }
    pub fn add_x_diff(&mut self, i: usize, x_diff: f64) {
        self.bins_for_even.add_x_diff(i, x_diff);
        self.bins_for_odd.add_x_diff(i, x_diff);
    }
}

fn log_polation(corr1: f64, corr2: f64, lag1: usize, lag2: usize, lag: usize) -> f64 {
    if (corr1 > corr2) && (corr2 > 0.0) {
        ((corr2.ln() - corr1.ln()) * ((lag - lag1) as f64)
            / ((lag2 - lag1) as f64)).exp()
    } else {
        0.0
    }
}

impl Covariances {
    fn new() -> Covariances {
        let short_lags = ShortLags::new();
        let long_lags = LongLags::new();
        Covariances { short_lags, long_lags }
    }
    fn add_x_diff(&mut self, i: usize, x_diff: f64) {
        self.short_lags.add_x_diff(x_diff);
        self.long_lags.bins_for_even.add_x_diff(i, x_diff);
        self.long_lags.bins_for_odd.add_x_diff(i, x_diff);
    }
    fn covariance(&self, n: usize, lag: usize) -> f64 {
        if lag >= n {
            0.0
        } else if lag <= SHORT_LAG_BUF_LENGTH {
            self.short_lags.cov_sums[lag - 1] / ((n - lag) as f64)
        } else if lag % 2 == 0 {
            let bins2 = &self.long_lags.bins_for_even;
            match Bins2::anchors_for(lag) {
                Anchors::Single(lag) =>
                    bins2.cov_sums[lag] / ((n / 2_usize.pow(lag as u32)) as f64),
                Anchors::Double(lag1, lag2) => {
                    let corr1 = bins2.cov_sums[lag1] / ((n / 2_usize.pow(lag as u32)) as f64);
                    let corr2 = bins2.cov_sums[lag2] / ((n / 2_usize.pow(lag as u32)) as f64);
                    log_polation(corr1, corr2, lag1, lag2, lag)
                }
            }
        } else {
            let bins3 = &self.long_lags.bins_for_odd;
            match Bins3::anchors_for(lag) {
                Anchors::Single(lag) =>
                    bins3.cov_sums[lag] / ((n / 3_usize.pow(lag as u32)) as f64),
                Anchors::Double(lag1, lag2) => {
                    let corr1 = bins3.cov_sums[lag1] / ((n / 3_usize.pow(lag as u32)) as f64);
                    let corr2 = bins3.cov_sums[lag2] / ((n / 3_usize.pow(lag as u32)) as f64);
                    log_polation(corr1, corr2, lag1, lag2, lag)
                }
            }
        }
    }
    pub(crate) fn add_covariance(&mut self, other: &Covariances) {
        self.short_lags.cov_sums
            .iter_mut()
            .zip(other.short_lags.cov_sums.iter())
            .for_each(|(sum, cov)| *sum += cov);
        self.long_lags.bins_for_even.cov_sums
            .iter_mut()
            .zip(other.long_lags.bins_for_even.cov_sums.iter())
            .for_each(|(sum, cov)| *sum += cov);
        self.long_lags.bins_for_odd.cov_sums
            .iter_mut()
            .zip(other.long_lags.bins_for_odd.cov_sums.iter())
            .for_each(|(sum, cov)| *sum += cov);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_thin_stats() {
        let mut stats = PolationStats::new();
        const N_ITERATIONS: usize = 1000;
        let mut a1: usize = 0;
        let mut a2: usize = 0;
        let mut a3: usize = 1;
        for i in 0..N_ITERATIONS {
            let a3_new = (a1 + a2 + a3) % 1024;
            a1 = a2;
            a2 = a3;
            a3 = a3_new;
            let x = (a3 as f64) + 1e4;
            stats.add(x);
            print!("{}\t{:.1}\t{:.1}\t{:.0}", i, x, stats.mean(), stats.variance());
            const LAG_MAX: usize = 24;
            for lag in 0..LAG_MAX {
                print!("\t{:.4}", stats.auto_corr(lag));
            }
            println!("\t{:.4}", stats.auto_corr_sum());
        }
    }
}