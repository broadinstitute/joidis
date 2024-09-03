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

struct Covariances {
    cov_sums: [f64; SHORT_LAG_BUF_LENGTH],
    x_diffs: VecDeque<f64>,
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
        self.n += 1;
        self.x_sum += x;
        let x_mean_previous = self.mean;
        self.mean = self.x_sum / (self.n as f64);
        let x_diff = x - self.mean;
        self.var_sum += (x - x_mean_previous) * x_diff;
        self.variance = self.var_sum / ((self.n - 1) as f64);
        self.covariances.add_x_diff(x_diff);
    }
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
    where
        I: Iterator<Item=&'a PolationStats>,
    {
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

impl Covariances {
    fn new() -> Covariances {
        let cov_sums = [0.0; SHORT_LAG_BUF_LENGTH];
        let x_diffs = VecDeque::with_capacity(SHORT_LAG_BUF_LENGTH);
        Covariances { cov_sums, x_diffs }
    }
    fn add_x_diff(&mut self, x_diff: f64) {
        for (j, x_diff_j) in self.x_diffs.iter().enumerate() {
            self.cov_sums[j] += x_diff * x_diff_j;
        }
        if self.x_diffs.len() == SHORT_LAG_BUF_LENGTH {
            self.x_diffs.pop_back();
        }
        self.x_diffs.push_front(x_diff);
    }
    fn covariance(&self, n: usize, lag: usize) -> f64 {
        if lag >= n {
            0.0
        } else if lag == 0 {
            1.0
        } else if lag <= SHORT_LAG_BUF_LENGTH {
            self.cov_sums[lag - 1] / ((n - lag) as f64)
        } else {
            let power = (lag / (SHORT_LAG_BUF_LENGTH + 1)) as i32;
            let lag = lag % (SHORT_LAG_BUF_LENGTH + 1);
            let cov_short = self.cov_sums[lag - 1] / ((n - lag) as f64);
            let cov_furthest =
                self.cov_sums[SHORT_LAG_BUF_LENGTH - 1] / ((n - SHORT_LAG_BUF_LENGTH) as f64);
            cov_short * cov_furthest.powi(power)
        }
    }
    pub(crate) fn add_covariance(&mut self, other: &Covariances) {
        self.cov_sums
            .iter_mut()
            .zip(other.cov_sums.iter())
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
            print!("{}\t{:.1}\t{:.1}\t{:.0}", i, x, stats.mean, stats.variance);
            const LAG_MAX: usize = 24;
            for lag in 0..LAG_MAX {
                print!("\t{:.4}", stats.auto_corr(lag));
            }
            println!("\t{:.4}", stats.auto_corr_sum());
        }
    }
}