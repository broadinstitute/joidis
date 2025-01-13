use rand::Rng;
use rand_distr::Normal;
use crate::params::Params;
use crate::sample::vars::Vars;
use rand_distr::Distribution;
use crate::data::GwasData;

pub(crate) struct GibbsSampler<R: Rng> {
    rng: R,
}

impl<R: Rng> GibbsSampler<R> {
    pub(crate) fn new(rng: R) -> GibbsSampler<R> { GibbsSampler { rng } }
    pub(crate) fn draw_e(&mut self, vars: &Vars, params: &Params, i_data_point: usize)
                         -> f64 {
        let n_traits = params.trait_names.len();
        let tau = params.tau;
        let inv_var_sum: f64 = 1.0 / tau.powi(2) + (0..n_traits).map(|i_trait| {
            (params.betas[i_trait] / params.sigmas[i_trait]).powi(2)
        }).sum::<f64>();
        let variance = 1.0 / inv_var_sum;
        let std_dev = variance.sqrt();
        let frac_sum = params.mu / tau.powi(2) + (0..n_traits).map(|i_trait| {
            params.betas[i_trait] * vars.ts[i_data_point][i_trait] / params.sigmas[i_trait].powi(2)
        }).sum::<f64>();
        let mean = variance * frac_sum;
        Normal::new(mean, std_dev).unwrap().sample(&mut self.rng)
    }

    pub(crate) fn draw_t(&mut self, data: &GwasData, vars: &Vars, params: &Params,
                         i_data_point: usize, i_trait: usize, pinned: bool) -> f64 {
        let mu_e = params.betas[i_trait] * vars.es[i_data_point];
        let var_e = params.sigmas[i_trait].powi(2);
        let mu_o = data.betas[i_data_point][i_trait];
        let se_o = data.ses[i_data_point][i_trait];
        if se_o > 0.0 && !pinned {
            let var_o = se_o.powi(2);
            let variance = 1.0 / (1.0 / var_e + 1.0 / var_o);
            let std_dev = variance.sqrt();
            let mean = variance * (mu_e / var_e + mu_o / var_o);
            Normal::new(mean, std_dev).unwrap().sample(&mut self.rng)
        } else {
            mu_o
        }
    }
}
