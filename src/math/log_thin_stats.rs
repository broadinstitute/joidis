pub(crate) struct LogThinStats {
    n: usize,
    sum: f64,
    m2: f64
}

impl LogThinStats {
    pub(crate) fn new() -> LogThinStats {
        let n: usize = 0;
        let sum: f64 = 0.0;
        let m2: f64 = 0.0;
        LogThinStats { n, sum, m2 }
    }
    pub(crate) fn add(&mut self, x: f64) {
        self.n += 1;
        self.sum += x;
    }
}