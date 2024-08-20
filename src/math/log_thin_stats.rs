pub(crate) struct LogThinStats {
    n: usize,
    x_sum: f64,
    var_sum: f64,
    mean: f64,
    variance: f64,
    betas: Betas
}

struct Betas {
    betas: Vec<f64>
}

impl LogThinStats {
    pub(crate) fn new() -> LogThinStats {
        let n: usize = 0;
        let x_sum: f64 = 0.0;
        let var_sum: f64 = 0.0;
        let mean: f64 = 0.0;
        let variance: f64 = 0.0;
        let betas = Betas::new();
        LogThinStats { n, x_sum, var_sum, mean, variance, betas }
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

impl Betas {
    pub(crate) fn new() -> Betas {
        Betas { betas: Vec::new() }
    }
    pub(crate) fn add(&mut self, beta: f64) {
        self.betas.push(beta);
    }
 }