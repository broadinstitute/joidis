use rand::Rng;
use crate::data::{GwasData, Meta};
use crate::sample::gibbs::GibbsSampler;
use crate::params::Params;
use crate::sample::var_stats::VarStats;
use crate::sample::vars::{VarIndex, Vars};

pub(crate) struct Sampler<R: Rng> {
    gibbs: GibbsSampler<R>,
    var_stats: VarStats,
}

pub(crate) trait ETracer {
    fn trace_e(&mut self, e: f64);
}

impl<R: Rng> Sampler<R> {
    pub(crate) fn new(meta: &Meta, rng: R) -> Sampler<R> {
        let gibbs = GibbsSampler::new(rng);
        let var_stats = VarStats::new(meta.clone());
        Sampler { gibbs, var_stats }
    }
    pub(crate) fn sample_n(&mut self, data: &GwasData, params: &Params, vars: &mut Vars,
                           n_steps: usize, e_tracer: &mut Option<Box<dyn ETracer>>,
                           t_pinned: bool) {
        for _ in 0..n_steps {
            self.sample_one(data, params, vars, e_tracer, t_pinned)
        }
    }
    pub(crate) fn sample_one(&mut self, data: &GwasData, params: &Params, vars: &mut Vars,
                             e_tracer: &mut Option<Box<dyn ETracer>>, t_pinned: bool) {
        for i_var in vars.indices() {
            match i_var {
                VarIndex::E { i_data_point } => {
                    let e = self.gibbs.draw_e(vars, params, i_data_point);
                    if let Some(e_tracer) = e_tracer {
                        e_tracer.trace_e(e);
                    }
                    vars.es[i_data_point] = e;
                }
                VarIndex::T { i_data_point, i_trait } => {
                    vars.ts[i_data_point][i_trait] =
                        self.gibbs.draw_t(data, vars, params, i_data_point, i_trait, t_pinned);
                }
            }
        }
        self.var_stats.add(vars);
    }
    pub(crate) fn var_stats(&self) -> &VarStats { &self.var_stats }
}
