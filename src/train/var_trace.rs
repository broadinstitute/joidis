use crate::math::wootz::WootzStats;
use crate::mcmc::metro::Draw;

pub(crate) struct VarTracer {
    stats: WootzStats,
}

const N_REFRESH_THRESHOLD: usize = 333;
const ATTEMPTS_REFRESH_THRESHOLD: usize = 14;

impl VarTracer {
    pub(crate) fn new(mean_estimate: f64, std_dev_estimate: f64) -> VarTracer {
        let x0 = mean_estimate - std_dev_estimate;
        let x1 = mean_estimate + std_dev_estimate;
        let stats = WootzStats::new(x0, x1);
        VarTracer { stats }
    }
    pub(crate) fn add(&mut self, draw: Draw) {
        self.stats.add(draw.x);
        if self.stats.n() > N_REFRESH_THRESHOLD &&
            draw.attempts_minus > ATTEMPTS_REFRESH_THRESHOLD &&
            draw.attempts_plus > ATTEMPTS_REFRESH_THRESHOLD {
            self.refresh();
        };
    }
    pub(crate) fn refresh(&mut self) {
        let mean = self.stats.mean();
        let std_dev = self.stats.variance().sqrt();
        let x0 = mean - std_dev;
        let x1 = mean + std_dev;
        self.stats = WootzStats::new(x0, x1);
    }
    pub(crate) fn std_dev(&self) -> f64 { self.stats.variance().sqrt() }
}