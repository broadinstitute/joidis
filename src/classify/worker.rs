use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender};
use log::{info, warn};
use rand::prelude::ThreadRng;
use rand::thread_rng;
use crate::classify::{Classification, MessageToCentral, MessageToWorker};
use crate::data::GwasData;
use crate::options::config::{ClassifyConfig, SharedConfig};
use crate::sample::vars::Vars;
use crate::params::Params;
use crate::sample::sampler::{Tracer, Sampler, NoOpTracer, defaults};
use crate::classify::exact::calculate_mu;
use crate::classify::tracer::ClassifyTracer;

pub(crate) fn classify_worker(data: &Arc<GwasData>, params: &Params, config: ClassifyConfig,
                              config_shared: SharedConfig, sender: Sender<MessageToCentral>,
                              receiver: Receiver<MessageToWorker>, i_thread: usize) {
    loop {
        let in_message = receiver.recv().unwrap();
        match in_message {
            MessageToWorker::DataPoint(i_data_point) => {
                let n_traits_total = data.meta.n_traits();
                let (data, is_col) = data.only_data_point(i_data_point);
                let trait_names = data.meta.trait_names.clone();
                if data.meta.n_traits() < n_traits_total {
                    let id = &data.meta.var_ids()[0];
                    let trait_list = trait_names.join(", ");
                    warn!("For {}, of the {} traits, we only have data for {} traits ({}).", id,
                        n_traits_total, data.meta.n_traits(), trait_list)
                }
                let params = params.reduce_to(trait_names, &is_col);
                let mut vars: Vec<Vars> =
                    (0..config.n_parallel())
                        .map(|_| Vars::initial_vars(&data, &params))
                        .collect();
                let rng = thread_rng();
                let meta = data.meta.clone();
                let n_chains = config.n_parallel();
                let mut sampler = Sampler::<ThreadRng>::new(&meta, rng, n_chains);
                let mut tracer =
                    match (&config.trace_ids, data.meta.var_ids.first()) {
                        (Some(trace_ids), Some(var_id))
                        if trace_ids.contains(var_id)
                        => {
                            let tracer =
                                ClassifyTracer::new(&meta, &config.out_file, var_id,
                                                    meta.n_endos(), n_chains);
                            Box::new(tracer) as Box<dyn Tracer>
                        }
                        _ => { Box::new(NoOpTracer::new()) }
                    };
                let n_steps_burn_in =
                    config_shared.n_steps_burn_in.unwrap_or(defaults::N_STEPS_BURN_IN);
                let var_ratio_burn_in =
                    config_shared.var_ratio_burn_in.unwrap_or(defaults::VAR_RATIO_BURN_IN);
                sampler.sample_n_ratio(&data, &params, &mut vars, n_steps_burn_in,
                                       var_ratio_burn_in, tracer.as_mut());
                sampler.reset_var_stats();
                let n_samples = config.n_samples.unwrap_or(defaults::N_SAMPLES);
                let var_ratio = config.var_ratio.unwrap_or(defaults::VAR_RATIO);
                sampler.sample_n_ratio(&data, &params, &mut vars, n_samples, var_ratio,
                                       tracer.as_mut());
                let sampled_classification =
                    sampler.var_stats().calculate_classification();
                let e_mean_calculated =
                    match calculate_mu(&params, &data.betas[0], &data.ses[0]) {
                        Ok(mu) => { mu }
                        Err(error) => {
                            info!("{}", error);
                            f64::NAN
                        }
                    };
                let classification =
                    Classification { i_data_point, sampled_classification, e_mean_calculated };
                sender.send(MessageToCentral { i_thread, classification }).unwrap();
            }
            MessageToWorker::Shutdown => {
                break;
            }
        }
    }
}