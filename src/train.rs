mod worker;
pub(crate) mod param_meta_stats;
mod initial_params;

use std::cmp;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender};
use std::thread::available_parallelism;
use crate::data::{GwasData, load_data};
use crate::error::Error;
use crate::options::action::Action;
use crate::options::cli::Flags;
use crate::options::config::{Config, SharedConfig};
use crate::report::Reporter;
use crate::train::initial_params::estimate_initial_params;
use crate::train::param_meta_stats::ParamMetaStats;
use crate::params::{Params, write_params_to_file};
use crate::sample::trace_file::ParamTraceFileWriter;
use crate::train::worker::train_worker;
use crate::util::display_option::DisplayOption;
use crate::util::threads::{InMessage, OutMessage, Threads, WorkerLauncher};

#[derive(Clone)]
pub(crate) enum MessageToWorker {
    TakeNSamples(usize),
    SetNewParams(Params),
    Shutdown,
}

impl OutMessage for MessageToWorker {
    const SHUTDOWN: Self = MessageToWorker::Shutdown;
}

pub(crate) struct MessageToCentral {
    i_thread: usize,
    params: Params,
}

impl MessageToCentral {
    pub(crate) fn new(i_thread: usize, params: Params) -> MessageToCentral {
        MessageToCentral { i_thread, params }
    }
}

impl InMessage for MessageToCentral {
    fn i_thread(&self) -> usize { self.i_thread }
}

#[derive(Clone)]
struct TrainWorkerLauncher {
    data: Arc<GwasData>,
    params: Params,
    config_shared: SharedConfig,
    t_pinned: bool,
}

impl WorkerLauncher<MessageToCentral, MessageToWorker> for TrainWorkerLauncher {
    fn launch(self, in_sender: Sender<MessageToCentral>, out_receiver: Receiver<MessageToWorker>,
              i_thread: usize) {
        let TrainWorkerLauncher {
            data, params, config_shared, t_pinned, ..
        } = self;
        train_worker(&data, params, in_sender, out_receiver, i_thread, &config_shared, t_pinned);
    }
}

impl TrainWorkerLauncher {
    fn new(data: Arc<GwasData>, params: Params, config_shared: SharedConfig, t_pinned: bool)
        -> TrainWorkerLauncher {
        TrainWorkerLauncher { data, params, config_shared, t_pinned }
    }
}

pub(crate) fn train_or_check(config: &Config, flags: Flags) -> Result<(), Error> {
    let data = load_data(config, Action::Train)?;
    println!("Loaded data for {} variants", data.meta.n_data_points());
    println!("{}", data);
    if flags.dry {
        println!("User picked dry run only, so doing nothing.")
    } else {
        train(data, config)?;
    }
    Ok(())
}

fn train(data: GwasData, config: &Config) -> Result<(), Error> {
    let data = Arc::new(data);
    let n_endos =  config.shared.n_endos;
    let n_traits = data.meta.n_traits();
    let mut params_trace_writer =
        if let Some(path) = &config.files.trace {
            let path = PathBuf::from(path);
            Some(ParamTraceFileWriter::new(path, n_endos, n_traits)?)
        } else {
            None
        };
    let n_threads = cmp::max(available_parallelism()?.get(), 3);
    println!("Launching {} workers. Configuration is: {} for iterations, {} for ratio.", n_threads,
             DisplayOption::new(config.shared.n_steps_burn_in),
             DisplayOption::new(config.shared.var_ratio_burn_in));
    let mut params = estimate_initial_params(&data, n_endos)?;
    println!("{}", params);
    let t_pinned = config.train.t_pinned.unwrap_or(false);
    let launcher =
        TrainWorkerLauncher::new(data, params.clone(), config.shared.clone(), t_pinned);
    let threads =
        Threads::<MessageToCentral, MessageToWorker>::new(launcher, n_threads);
    println!("Workers launched and burned in.");
    let n_samples: usize = config.train.n_samples_per_iteration;
    let mut reporter = Reporter::new();
    let mut i_round: usize = 0;
    let mut i_iteration: usize = 0;
    loop {
        let params0 = create_param_estimates(&threads, n_samples)?;
        let params1 = create_param_estimates(&threads, n_samples)?;
        let mut param_meta_stats =
            ParamMetaStats::new(n_threads, n_endos, params.trait_names.clone(),
                                &params0, &params1);
        let mut reached_precision = false;
        loop {
            i_iteration += 1;
            let params_new = create_param_estimates(&threads, n_samples)?;
            param_meta_stats.add(&params_new);
            let summary = param_meta_stats.summary()?;
            if i_iteration >= config.train.n_iterations_per_round {
                params = summary.params.clone();
                if let Some(params_trace_writer) = &mut params_trace_writer {
                    params_trace_writer.write(&params)?;
                }
                if i_round >= config.train.n_rounds {
                    println!("Done!");
                    reached_precision = true;
                } else {
                    i_round += 1;
                    println!("Setting new parameters for round {} after {} iterations", i_round,
                             i_iteration);
                    i_iteration = 0;
                    for sender in threads.out_senders.iter() {
                        sender.send(MessageToWorker::SetNewParams(params.clone()))?;
                    }
                }
                reporter.report(&summary, i_round, i_iteration, n_samples);
                reporter.reset_round_timer();
                break;
            }
        }
        if reached_precision {
            break;
        }
    };
    if config.train.normalize_mu_to_one {
        params = params.normalized_with_mu_one()
    }
    write_params_to_file(&params, config.files.params.as_str())?;
    Ok(())
}

fn create_param_estimates(threads: &Threads<MessageToCentral, MessageToWorker>, n_samples: usize)
                          -> Result<Vec<Params>, Error> {
    threads.broadcast(MessageToWorker::TakeNSamples(n_samples))?;
    let responses = threads.responses_from_all()?;
    let params: Vec<Params> =
        responses.into_iter().map(|response| response.params).collect();
    Ok(params)
}
