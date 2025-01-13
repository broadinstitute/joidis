use std::fs::read_to_string;
use serde::{Deserialize, Serialize};
use crate::error::{Error, for_file};
use crate::params::ParamsOverride;
use crate::data::gwas::GwasCols;

#[derive(Deserialize, Serialize)]
pub(crate) struct Config {
    pub(crate) files: FilesConfig,
    pub(crate) gwas: Vec<GwasConfig>,
    pub(crate) train: TrainConfig,
    pub(crate) classify: ClassifyConfig,
}

#[derive(Deserialize, Serialize, Clone)]
pub(crate) struct GwasConfig {
    pub(crate) name: String,
    pub(crate) file: String,
    pub(crate) cols: Option<GwasCols>
}

#[derive(Deserialize, Serialize)]
pub(crate) struct FilesConfig {
    pub(crate) trace: Option<String>,
    pub(crate) params: String
}

#[derive(Deserialize, Serialize, Clone)]
pub(crate) struct TrainConfig {
    pub(crate) ids_file: String,
    pub(crate) n_steps_burn_in: usize,
    pub(crate) n_samples_per_iteration: usize,
    pub(crate) n_iterations_per_round: usize,
    pub(crate) n_rounds: usize,
    pub(crate) normalize_mu_to_one: bool,
    pub(crate) params_trace_file: Option<String>,
    pub(crate) t_pinned: Option<bool>,
}

#[derive(Deserialize, Serialize, Clone)]
pub(crate) struct ClassifyConfig {
    pub(crate) params_override: Option<ParamsOverride>,
    pub(crate) n_steps_burn_in: usize,
    pub(crate) n_samples: usize,
    pub(crate) out_file: String,
    pub(crate) trace_ids: Option<Vec<String>>,
    pub(crate) t_pinned: Option<bool>,
}

pub(crate) fn load_config(file: &str) -> Result<Config, Error> {
    let string = for_file(file, read_to_string(file))?;
    let config: Config = toml::from_str(&string)?;
    Ok(config)
}
