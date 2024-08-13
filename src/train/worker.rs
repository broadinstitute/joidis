use std::slice;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender};
use rand::prelude::ThreadRng;
use rand::thread_rng;
use crate::data::GwasData;
use crate::options::config::SharedConfig;
use crate::train::{MessageToCentral, MessageToWorker};
use crate::params::Params;
use crate::sample::sampler::{defaults, NoOpTracer, Sampler};
use crate::sample::vars::Vars;

pub(crate) fn train_worker(data: &Arc<GwasData>, mut params: Params,
                           sender: Sender<MessageToCentral>, receiver: Receiver<MessageToWorker>,
                           i_thread: usize, config_shared: &SharedConfig) {
    let mut vars = Vars::initial_vars(data, &params);
    let rng = thread_rng();
    let meta = data.meta.clone();
    let n_chains: usize = 1;
    let mut sampler = Sampler::<ThreadRng>::new(&meta, rng, n_chains);
    let n_steps_burn_in = config_shared.n_steps_burn_in.unwrap_or(defaults::N_STEPS_BURN_IN);
    let var_ratio_burn_in =
        config_shared.var_ratio_burn_in.unwrap_or(defaults::VAR_RATIO_BURN_IN);
    let mut tracer = NoOpTracer::new();
    sampler.sample_n_ratio(data, &params, slice::from_mut(&mut vars), n_steps_burn_in,
                           var_ratio_burn_in, &mut tracer);
    loop {
        let in_message = receiver.recv().unwrap();
        match in_message {
            MessageToWorker::TakeNSamples(n_samples) => {
                sampler.sample_n(data, &params, slice::from_mut(&mut vars), n_samples,
                                 &mut tracer);
                let params_new = sampler.var_stats().compute_new_params();
                sender
                    .send(MessageToCentral::new(i_thread, params_new))
                    .unwrap();
            }
            MessageToWorker::SetNewParams(params_new) => {
                params = params_new;
                sampler.sample_n_ratio(data, &params, slice::from_mut(&mut vars), n_steps_burn_in,
                                       var_ratio_burn_in, &mut tracer);
            }
            MessageToWorker::Shutdown => {
                break;
            }
        }
    }
}
