use crate::data::Meta;
use crate::math::matrix::Matrix;
use crate::params::Params;
use crate::sample::vars::Vars;

pub(crate) struct VarStats {
    meta: Meta,
    n: usize,
    e_sums: Matrix,
    e2_sums: Matrix,
    // e_t_sums: Matrix,
    t_sums: Matrix,
    t2_sums: Matrix,
}

pub(crate) struct SampledClassification {
    pub(crate) e_mean: Vec<f64>,
    pub(crate) e_std: Vec<f64>,
    pub(crate) t_means: Vec<f64>,
}

impl VarStats {
    pub(crate) fn new(meta: Meta) -> VarStats {
        let n: usize = 0;
        let n_data_points = meta.n_data_points();
        let n_endos = meta.n_endos();
        let n_traits = meta.n_traits();
        let e_sums: Matrix = Matrix::fill(n_data_points, n_endos, |_, _| 0.0);
        let e2_sums: Matrix = Matrix::fill(n_data_points, n_endos, |_, _| 0.0);
        // let e_t_sums: Matrix = Matrix::fill(n_data_points, n_traits, |_, _| 0.0);
        let t_sums: Matrix = Matrix::fill(n_data_points, n_traits, |_, _| 0.0);
        let t2_sums: Matrix = Matrix::fill(n_data_points, n_traits, |_, _| 0.0);
        VarStats { meta, n, e_sums, e2_sums, /* e_t_sums, */ t_sums, t2_sums }
    }
    pub(crate) fn add(&mut self, vars: &Vars) {
        self.n += 1;
        let n_data_points = self.meta.n_data_points();
        let n_traits = self.meta.n_traits();
        let n_endos = vars.meta.n_endos();
        for i_data_point in 0..n_data_points {
            for i_endo in 0..n_endos {
                let e_j_k = vars.es[i_data_point][i_endo];
                self.e_sums[i_data_point][i_endo] += e_j_k;
                self.e2_sums[i_data_point][i_endo] += e_j_k.powi(2);
            }
            for i_trait in 0..n_traits {
                let t_j_i = vars.ts[i_data_point][i_trait];
                // self.e_t_sums[i_data_point][i_trait] += e_j * t_j_i;
                self.t_sums[i_data_point][i_trait] += t_j_i;
                self.t2_sums[i_data_point][i_trait] += t_j_i.powi(2);
            }
        }
    }
    pub(crate) fn sum(stats_list: &[VarStats]) -> VarStats {
        let meta = stats_list[0].meta.clone();
        let n = stats_list.iter().map(|stats| stats.n).sum();
        let n_data_points = meta.n_data_points();
        let n_endos = meta.n_endos();
        let n_traits = meta.n_traits();
        let e_sums = Matrix::fill(n_data_points, n_endos, |i, j| {
            stats_list.iter().map(|stats| stats.e_sums[i][j]).sum()
        });
        let e2_sums = Matrix::fill(n_data_points, n_endos, |i, j| {
            stats_list.iter().map(|stats| stats.e2_sums[i][j]).sum()
        });
        // let e_t_sums = Matrix::fill(n_data_points, n_traits, |i, j| {
        //     stats_list.iter().map(|stats| stats.e_t_sums[i][j]).sum()
        // });
        let t_sums = Matrix::fill(n_data_points, n_traits, |i, j| {
            stats_list.iter().map(|stats| stats.t_sums[i][j]).sum()
        });
        let t2_sums = Matrix::fill(n_data_points, n_traits, |i, j| {
            stats_list.iter().map(|stats| stats.t2_sums[i][j]).sum()
        });
        VarStats { meta, n, e_sums, e2_sums, /* e_t_sums, */ t_sums, t2_sums }
    }
    pub(crate) fn compute_new_params(&self) -> Params {
        // let meta = &self.meta;
        // let n_f = self.n as f64;
        // let n_data_points = meta.n_data_points();
        // let n_data_points_f = n_data_points as f64;
        // let n_traits = meta.n_traits();
        // let mut sum_for_mu: f64 = 0.0;
        // for j in 0..n_data_points {
        //     let mean_e_j = self.e_sums[j] / n_f;
        //     sum_for_mu += mean_e_j;
        // }
        // let mu = sum_for_mu / n_data_points_f;
        // let mut sum_for_tau: f64 = 0.0;
        // for j in 0..n_data_points {
        //     let mean_e2_j = self.e2_sums[j] / n_f;
        //     let mean_e_j = self.e_sums[j] / n_f;
        //     sum_for_tau += mean_e2_j - 2.0 * mu * mean_e_j + n_data_points_f * mu.powi(2);
        // }
        // let tau = (sum_for_tau / n_data_points_f).sqrt();
        // let mut betas: Vec<f64> = Vec::with_capacity(n_traits);
        // for i in 0..n_traits {
        //     let mut mean_e_t_sum: f64 = 0.0;
        //     let mut mean_e2_sum: f64 = 0.0;
        //     for j in 0..n_data_points {
        //         mean_e_t_sum += self.e_t_sums[j][i] / n_f;
        //         mean_e2_sum += self.e2_sums[j] / n_f;
        //     }
        //     betas.push(mean_e_t_sum / mean_e2_sum);
        // }
        // let mut sigmas: Vec<f64> = Vec::with_capacity(n_traits);
        // for (i, beta) in betas.iter().enumerate() {
        //     let mut sum_for_sigma: f64 = 0.0;
        //     for j in 0..n_data_points {
        //         let mean_t2_j_i = self.t2_sums[j][i] / n_f;
        //         let mean_e_t_j_i = self.e_t_sums[j][i] / n_f;
        //         let mean_e2_j_i = self.e2_sums[j] / n_f;
        //         sum_for_sigma +=
        //             mean_t2_j_i - 2.0 * betas[i] * mean_e_t_j_i + beta.powi(2) * mean_e2_j_i
        //     }
        //     let sigma = (sum_for_sigma / n_data_points_f).sqrt();
        //     sigmas.push(sigma)
        // }
        // let trait_names = meta.trait_names.clone();
        // Params { trait_names, mu, tau, betas, sigmas }
        todo!()
    }
    pub(crate) fn calculate_classification(&self) -> SampledClassification {
        let meta = &self.meta;
        let n_data_points = meta.n_data_points();
        let n_endos = meta.n_endos();
        let n_traits = meta.n_traits();
        let denom = (self.n * n_data_points) as f64;
        let mut e_mean: Vec<f64> = vec![0.0; n_endos];
        let mut e2_mean: Vec<f64> = vec![0.0; n_endos];
        let mut t_means: Vec<f64> = vec![0.0; n_traits];
        for i_data_point in 0..n_data_points {
            for i_endo in 0..n_endos {
                e_mean[i_endo] += self.e_sums[i_data_point][i_endo] / denom;
                e2_mean[i_endo] += self.e2_sums[i_data_point][i_endo] / denom;
            }
            for (i, t_mean) in t_means.iter_mut().enumerate() {
                *t_mean += self.t_sums[i_data_point][i] / denom;
            }
        }
        let e_std =
            (0..n_endos).map(|i_endo|
            (e2_mean[i_endo] - e_mean[i_endo].powi(2)).sqrt()).collect::<Vec<f64>>();
        SampledClassification { e_mean, e_std, t_means }
    }
    pub(crate) fn calculate_convergences(var_stats_list: &[VarStats])
        -> impl Iterator<Item=f64> + '_{
        let meta = &var_stats_list[0].meta;
        let n_data_points = meta.n_data_points();
        let n_endos = meta.n_endos();
        let n_traits = meta.n_traits();
        let e_ratio_iter = (0..n_endos).map(move |i_endo| {
            (0..n_data_points).map(|i_data_point| {
                let e_var_mean: f64 =
                    var_stats_list.iter().map(|stats| {
                        let n = stats.n as f64;
                        let e_mean = stats.e_sums[i_data_point][i_endo] / n;
                        let e2_mean = stats.e2_sums[i_data_point][i_endo] / n;
                        e2_mean - e_mean.powi(2)
                    }).sum::<f64>() / (var_stats_list.len() as f64);
                let e_mean_mean: f64 =
                    var_stats_list.iter().map(|stats| {
                        let n = stats.n as f64;
                        stats.e_sums[i_data_point][i_endo] / n
                    }).sum::<f64>() / (var_stats_list.len() as f64);
                let e_mean_var: f64 =
                    var_stats_list.iter().map(|stats| {
                        let n = stats.n as f64;
                        let e_mean: f64 = stats.e_sums[i_data_point][i_endo] / n;
                        (e_mean - e_mean_mean).powi(2)
                    }).sum::<f64>() / (var_stats_list.len() as f64);
                e_mean_var / e_var_mean
            }).sum::<f64>() / (n_data_points as f64)
        });
        let t_ratio_iter = (0..n_traits).map(move |i_trait| {
            (0..n_data_points).map(|i_data_point| {
                let t_var_mean: f64 =
                    var_stats_list.iter().map(|stats| {
                        let n = stats.n as f64;
                        let t_mean = stats.t_sums[i_data_point][i_trait] / n;
                        let t2_mean = stats.t2_sums[i_data_point][i_trait] / n;
                        t2_mean - t_mean.powi(2)
                    }).sum::<f64>() / (var_stats_list.len() as f64);
                let t_mean_mean: f64 =
                    var_stats_list.iter().map(|stats| {
                        let n = stats.n as f64;
                        stats.t_sums[i_data_point][i_trait] / n
                    }).sum::<f64>() / (var_stats_list.len() as f64);
                let t_mean_var: f64 =
                    var_stats_list.iter().map(|stats| {
                        let n = stats.n as f64;
                        let t_mean = stats.t_sums[i_data_point][i_trait] / n;
                        (t_mean - t_mean_mean).powi(2)
                    }).sum::<f64>() / (var_stats_list.len() as f64);
                t_mean_var / t_var_mean
            }).sum::<f64>() / (n_data_points as f64)
        });
        e_ratio_iter.chain(t_ratio_iter)
    }
    pub(crate) fn reset(&mut self) {
        self.n = 0;
        self.e_sums.elements.fill(0.0);
        self.e2_sums.elements.fill(0.0);
        self.t_sums.elements.fill(0.0);
        self.t2_sums.elements.fill(0.0);
    }
}