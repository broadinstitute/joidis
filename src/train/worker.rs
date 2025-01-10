use std::slice;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender};
use rand::prelude::ThreadRng;
use rand::thread_rng;
use crate::data::GwasData;
use crate::options::config::SharedConfig;
use crate::train::{MessageToCentral, MessageToWorker};
use crate::params::Params;
use crate::sample::sampler::{NoOpTracer, Sampler, StopConditions};
use crate::sample::vars::Vars;

pub(crate) fn train_worker(data: &Arc<GwasData>, mut params: Params,
                           sender: Sender<MessageToCentral>, receiver: Receiver<MessageToWorker>,
                           i_thread: usize, config_shared: &SharedConfig) {
    let mut vars = Vars::initial_vars(data, &params);
    let rng = thread_rng();
    let meta = data.meta.clone();
    let n_chains: usize = 1;
    let mut sampler = Sampler::<ThreadRng>::new(&meta, rng, n_chains);
    let mut tracer = NoOpTracer::new();
    let stop_conditions = StopConditions::for_burn_in(config_shared);
    let t_pinned = config_shared.t_pinned.unwrap_or(false);
    sampler.sample_conditional(data, &params, slice::from_mut(&mut vars), &stop_conditions,
                               &mut tracer, t_pinned);
    loop {
        let in_message = receiver.recv().unwrap();
        match in_message {
            MessageToWorker::TakeNSamples(n_samples) => {
                sampler.sample_n(data, &params, slice::from_mut(&mut vars), n_samples,
                                 &mut tracer, t_pinned);
                let params_new = sampler.var_stats().compute_new_params();
                sender
                    .send(MessageToCentral::new(i_thread, params_new))
                    .unwrap();
            }
            MessageToWorker::SetNewParams(params_new) => {
                params = params_new;
                sampler.sample_conditional(data, &params, slice::from_mut(&mut vars),
                                           &stop_conditions, &mut tracer, t_pinned);
            }
            MessageToWorker::Shutdown => {
                break;
            }
        }
    }
}
