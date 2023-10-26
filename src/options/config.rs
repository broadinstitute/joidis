use std::fs::read_to_string;
use serde::Deserialize;
use crate::error::Error;

#[derive(Deserialize)]
pub(crate) struct Config {
    pub(crate) gwas: Vec<GwasConfig>,
    pub(crate) train: TrainConfig,
}

#[derive(Deserialize, Clone)]
pub(crate) struct TrainConfig {
    pub(crate) ids_file: String,
    pub(crate) n_steps_per_sample: usize,
    pub(crate) n_steps_burn_in: usize,
    pub(crate) n_samples_per_round: usize,
    pub(crate) precision: f64
}

#[derive(Deserialize)]
pub(crate) struct GwasConfig {
    pub(crate) name: String,
    pub(crate) file: String
}

pub(crate) fn load_config(file: &str) -> Result<Config, Error> {
    let string = read_to_string(file)?;
    let config: Config = toml::from_str(&string)?;
    Ok(config)
}
