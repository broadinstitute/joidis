use crate::data::{GwasData, Meta};
use crate::params::Params;
use crate::sample::gibbs::GibbsSampler;
use crate::sample::var_stats::VarStats;
use crate::sample::vars::{VarIndex, Vars};
use rand::Rng;
use crate::options::config::{ClassifyConfig, SharedConfig};

pub(crate) mod defaults {
    pub(crate) const N_STEPS_BURN_IN: usize = 1000;
    pub(crate) const N_STEPS_EFFECTIVE_BURN_IN: usize = 100;
    pub(crate) const VAR_RATIO_BURN_IN: f64 = 0.1;
    pub(crate) const N_SAMPLES: usize = 2000;
    pub(crate) const N_SAMPLES_EFFECTIVE: usize = 20;
    pub(crate) const VAR_RATIO: f64 = 0.05;
}

pub(crate) struct StopConditions {
    pub(crate) n_samples: usize,
    pub(crate) n_samples_effective: usize,
    pub(crate) var_ratio: f64,
}

impl StopConditions {
    pub(crate) fn new(n_samples: usize, n_samples_effective: usize, var_ratio: f64) -> StopConditions {
        StopConditions { n_samples, n_samples_effective, var_ratio }
    }
    pub(crate) fn for_burn_in(shared_config: &SharedConfig) -> StopConditions {
        let n_steps = shared_config.n_steps_burn_in.unwrap_or(defaults::N_STEPS_BURN_IN);
        let n_steps_effective =
            shared_config.n_steps_effective_burn_in.unwrap_or(defaults::N_STEPS_EFFECTIVE_BURN_IN);
        let var_ratio =
            shared_config.var_ratio_burn_in.unwrap_or(defaults::VAR_RATIO_BURN_IN);
        StopConditions::new(n_steps, n_steps_effective, var_ratio)
    }
    pub(crate) fn for_classification(classify_config: &ClassifyConfig) -> StopConditions {
        let n_samples = classify_config.n_samples.unwrap_or(defaults::N_SAMPLES);
        let n_samples_effective =
            classify_config.n_samples_effective.unwrap_or(defaults::N_SAMPLES_EFFECTIVE);
        let var_ratio = classify_config.var_ratio.unwrap_or(defaults::VAR_RATIO);
        StopConditions::new(n_samples, n_samples_effective, var_ratio)
    }
}

pub(crate) struct Sampler<R: Rng> {
    gibbs: GibbsSampler<R>,
    var_stats_list: Vec<VarStats>,
}

pub(crate) trait Tracer {
    fn trace_e(&mut self, i_endo: usize, e: f64, i_chain: usize);
    fn trace_t(&mut self, i_trait: usize, t: f64, i_chain: usize);
    fn write_values(&mut self);
    fn trace_convergence(&mut self, var_stats_list: &[VarStats]);
    fn trace_n_effective(&mut self, var_stats_list: &[VarStats]);
}

pub(crate) struct NoOpTracer;

impl NoOpTracer {
    pub(crate) fn new() -> NoOpTracer { NoOpTracer }
}

impl Tracer for NoOpTracer {
    fn trace_e(&mut self, _i_endo: usize, _e: f64, _i_chain: usize) {}
    fn trace_t(&mut self, _i_trait: usize, _t: f64, _i_chain: usize) {}

    fn write_values(&mut self) {}

    fn trace_convergence(&mut self, _var_stats_list: &[VarStats]) {}

    fn trace_n_effective(&mut self, _var_stats_list: &[VarStats]) {}
}

impl<R: Rng> Sampler<R> {
    pub(crate) fn new(meta: &Meta, rng: R, n_chains: usize) -> Sampler<R> {
        let gibbs = GibbsSampler::new(rng);
        let var_stats_list =
            (0..n_chains).map(|_| VarStats::new(meta.clone())).collect();
        Sampler { gibbs, var_stats_list }
    }
    pub(crate) fn sample_conditional(&mut self, data: &GwasData, params: &Params, vars: &mut [Vars],
                                     stop_conditions: &StopConditions, tracer: &mut dyn Tracer) {
        self.sample_n(data, params, vars, stop_conditions.n_samples, tracer);
        let mut n_steps2 = 10;
        loop {
            let max_convergence =
                VarStats::calculate_convergences(&self.var_stats_list)
                    .fold(f64::NEG_INFINITY, f64::max);
            let n_effective =
                self.var_stats_list.iter().map(VarStats::effective_sample_size)
                    .fold(f64::INFINITY, f64::min);
            if max_convergence < stop_conditions.var_ratio
                && n_effective >=stop_conditions.n_samples_effective as f64 {
                break;
            } else {
                n_steps2 += n_steps2 / 10 + 1;
            }
            self.sample_n(data, params, vars, n_steps2, tracer);
        }
    }
    pub(crate) fn sample_n(&mut self, data: &GwasData, params: &Params, vars: &mut [Vars],
                           n_steps: usize, tracer: &mut dyn Tracer) {
        for _ in 0..n_steps {
            self.sample_one(data, params, vars, tracer)
        }
    }
    pub(crate) fn sample_one(&mut self, data: &GwasData, params: &Params, vars: &mut [Vars],
                             tracer: &mut dyn Tracer) {
        for (i_chain, vars) in vars.iter_mut().enumerate() {
            for i_var in vars.indices() {
                match i_var {
                    VarIndex::E { i_data_point, i_endo } => {
                        let e = self.gibbs.draw_e(vars, params, i_data_point, i_endo);
                        tracer.trace_e(i_endo, e, i_chain);
                        vars.es[i_data_point][i_endo] = e;
                    }
                    VarIndex::T { i_data_point, i_trait } => {
                        let t = self.gibbs.draw_t(data, vars, params, i_data_point, i_trait);
                        tracer.trace_t(i_trait, t, i_chain);
                        vars.ts[i_data_point][i_trait] = t;
                    }
                }
            }
            self.var_stats_list[i_chain].add(vars);
        }
        tracer.write_values();
        tracer.trace_convergence(&self.var_stats_list);
        tracer.trace_n_effective(&self.var_stats_list);
    }
    pub(crate) fn var_stats(&self) -> VarStats { VarStats::sum(&self.var_stats_list) }
    pub(crate) fn reset_var_stats(&mut self) {
        self.var_stats_list.iter_mut().for_each(|stats| stats.reset());
    }
}
