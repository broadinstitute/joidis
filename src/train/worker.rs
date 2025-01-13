use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender};
use rand::prelude::ThreadRng;
use rand::thread_rng;
use crate::data::LoadedData;
use crate::options::config::TrainConfig;
use crate::train::{MessageToCentral, MessageToWorker};
use crate::params::Params;
use crate::sample::sampler::Sampler;
use crate::sample::vars::Vars;

pub(crate) fn train_worker(data: &Arc<LoadedData>, mut params: Params,
                           sender: Sender<MessageToCentral>, receiver: Receiver<MessageToWorker>,
                           i_thread: usize, config: &TrainConfig) {
    let mut vars = Vars::initial_vars(&data.gwas_data, &params);
    let rng = thread_rng();
    let meta = data.gwas_data.meta.clone();
    let mut sampler = Sampler::<ThreadRng>::new(&meta, rng);
    let t_pinned = config.t_pinned.unwrap_or(false);
    sampler.sample_n(&data.gwas_data, &params, &mut vars, config.n_steps_burn_in, &mut None,
                     t_pinned);
    loop {
        let in_message = receiver.recv().unwrap();
        match in_message {
            MessageToWorker::TakeNSamples(n_samples) => {
                sampler.sample_n(&data.gwas_data, &params, &mut vars, n_samples, &mut None,
                                 t_pinned);
                let params_new = sampler.var_stats().compute_new_params(&data.weights);
                sender
                    .send(MessageToCentral::new(i_thread, params_new))
                    .unwrap();
            }
            MessageToWorker::SetNewParams(params_new) => {
                params = params_new;
                sampler.sample_n(&data.gwas_data, &params, &mut vars, config.n_steps_burn_in,
                                 &mut None, t_pinned);
            }
            MessageToWorker::Shutdown => {
                break;
            }
        }
    }
}
