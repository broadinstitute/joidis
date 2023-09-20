mod model;
mod param_stats;
mod sampler;
mod params;
mod vars;

use crate::data::{load_training_data, TrainData};
use crate::error::Error;
use crate::options::config::Config;
use crate::train::model::TrainModel;
use crate::train::param_stats::ParamDiffStats;
use crate::train::params::Params;
use crate::train::sampler::Sampler;

pub(crate) fn train_or_check(config: &Config, dry: bool) -> Result<(), Error> {
    let data = load_training_data(config)?;
    println!("Loaded data for {} variants", data.beta_se_lists.len());
    println!("{}", data);
    if !dry {
        train(data)?;
    }
    Ok(())
}

fn train(data: TrainData) -> Result<Params, Error> {
    let model = TrainModel::new(data);
    let mut params = model.initial_params();
    let mut vars = model.initial_vars(&params);
    let sampler = Sampler::new();
    loop {
        let mut stats = ParamDiffStats::new();
        let stats =
            loop {
                sampler.sample(&model, &params, &mut vars);
                let sample = model.evaluate_params(&params, &vars);
                stats.add_diffs(sample);
                if stats.ready_for_param_estimate() {
                    break stats;
                }
            };
        let estimate = stats.estimate_params();
        params = estimate.params;
        if estimate.is_done { break }
    }
    Ok(params)
}