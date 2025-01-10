use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;
use crate::data::gwas::GwasCols;
use crate::error::{Error, for_file};
use crate::math::matrix::Matrix;
use crate::options::cli::ImportPhenetOptions;
use crate::options::config::{ClassifyConfig, Config, FilesConfig, GwasConfig, SharedConfig, TrainConfig, write_config};
use crate::params::{Params, ParamsOverride};

mod defaults {
    pub(crate) mod shared {
        pub(crate) const N_STEPS_BURN_IN: usize = 10000;
        pub(crate) const N_STEPS_EFFECTIVE_BURN_IN: usize = 100;
        pub(crate) const VAR_RATIO_BURN_IN: f64 = 0.1;
    }
    pub(crate) mod train {
        pub(crate) const N_SAMPLES_PER_ITERATION: usize = 100;
        pub(crate) const N_ITERATIONS_PER_ROUND: usize = 1000;
        pub(crate) const N_ROUNDS: usize = 10000;
    }

    pub(crate) mod classify {
        pub(crate) const N_SAMPLES: usize = 100_000;
        pub(crate) const N_SAMPLES_EFFECTIVE: usize = 1000;
        pub(crate) const VAR_RATIO: f64 = 0.05;
        pub(crate) const N_PARALLEL: Option<usize> = None;
    }
}

mod keys {
    pub(crate) const VAR_ID_FILE: &str = "var_id_file";
    pub(crate) const CONFIG_FILE: &str = "config_file";
    pub(crate) const OUTPUT_FILE: &str = "output_file";
    pub(crate) const DECLARE: &str = "declare";
    pub(crate) const TRAIT: &str = "trait";
    pub(crate) const ENDO: &str = "endo";
    pub(crate) const FILE: &str = "file";
    pub(crate) const ID_COL: &str = "id_col";
    pub(crate) const EFFECT_COL: &str = "effect_col";
    pub(crate) const SE_COL: &str = "se_col";
    pub(crate) const BETA: &str = "beta";
    pub(crate) const VAR: &str = "var";
    pub(crate) const MEAN: &str = "mean";
}

struct PhenetOpts {
    var_id_file: String,
    config_files: Vec<String>,
    output_file: Option<String>,
}

struct ConfigBuilder {
    trait_names: Vec<String>,
    endo_names: Vec<String>,
    files: BTreeMap<String, String>,
    id_cols: BTreeMap<String, String>,
    effect_cols: BTreeMap<String, String>,
    se_cols: BTreeMap<String, String>,
    betas: BTreeMap<String, f64>,
    vars: BTreeMap<String, f64>,
    means: BTreeMap<String, f64>,
}

impl ConfigBuilder {
    fn new() -> ConfigBuilder {
        let trait_names: Vec<String> = Vec::new();
        let endo_names: Vec<String> = Vec::new();
        let files: BTreeMap<String, String> = BTreeMap::new();
        let id_cols: BTreeMap<String, String> = BTreeMap::new();
        let effect_cols: BTreeMap<String, String> = BTreeMap::new();
        let se_cols: BTreeMap<String, String> = BTreeMap::new();
        let betas: BTreeMap<String, f64> = BTreeMap::new();
        let vars: BTreeMap<String, f64> = BTreeMap::new();
        let means: BTreeMap<String, f64> = BTreeMap::new();
        ConfigBuilder {
            trait_names,
            endo_names,
            files,
            id_cols,
            effect_cols,
            se_cols,
            betas,
            vars,
            means,
        }
    }
    fn read_phenet_config(&mut self, file: &str) -> Result<(), Error> {
        let file = for_file(file, File::open(file))?;
        let reader = BufReader::new(file);
        self.parse_phenet_config(reader)?;
        Ok(())
    }
    fn read_phenet_config_optional(&mut self, file: &str) -> Result<(), Error> {
        match for_file(file, File::open(file)) {
            Ok(file) => {
                let reader = BufReader::new(file);
                self.parse_phenet_config(reader)?;
            }
            Err(error) => {
                println!("{}", error);
                println!("Since this file was optional, we proceed.");
            }
        }
        Ok(())
    }
    fn parse_phenet_config(&mut self, reader: BufReader<File>) -> Result<(), Error> {
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if !line.is_empty() {
                let mut parts_iter = line.split(char::is_whitespace);
                match (parts_iter.next(), parts_iter.next(), parts_iter.next(), parts_iter.next()) {
                    (Some(part1), Some(part2), Some(part3), None) => {
                        match (part1, part2, part3) {
                            (trait_name, keys::DECLARE, keys::TRAIT) => {
                                let trait_name = trait_name.to_string();
                                if !self.trait_names.contains(&trait_name) {
                                    self.trait_names.push(trait_name.to_string());
                                }
                            }
                            (endo_name, keys::DECLARE, keys::ENDO) => {
                                self.endo_names.push(endo_name.to_string());
                            }
                            (trait_name, keys::FILE, file) => {
                                self.files.insert(trait_name.to_string(),
                                                  file.to_string());
                            }
                            (trait_name, keys::ID_COL, id_col) => {
                                self.id_cols.insert(trait_name.to_string(),
                                                    id_col.to_string());
                            }
                            (trait_name, keys::EFFECT_COL, effect_col) => {
                                self.effect_cols.insert(trait_name.to_string(),
                                                        effect_col.to_string());
                            }
                            (trait_name, keys::SE_COL, se_col) => {
                                self.se_cols.insert(trait_name.to_string(),
                                                    se_col.to_string());
                            }
                            (trait_name, keys::BETA, beta) => {
                                self.betas.insert(trait_name.to_string(),
                                                  beta.parse::<f64>()?);
                            }
                            (trait_name, keys::VAR, var) => {
                                self.vars.insert(trait_name.to_string(),
                                                 var.parse::<f64>()?);
                            }
                            (trait_name, keys::MEAN, mean) => {
                                self.means.insert(trait_name.to_string(),
                                                  mean.parse::<f64>()?);
                            }
                            _ => {
                                println!("Ignoring line '{}'.", line)
                            }
                        }
                    }
                    _ => {
                        Err(Error::from(format!("Cannot parse line: '{}'.", line)))?;
                    }
                }
            }
        }
        Ok(())
    }
    fn report(&self) {
        println!("Trait names: {}", self.trait_names.join(", "));
        println!("Endophenotype names: {}", self.endo_names.join(", "));
        let file_default = "<no file>".to_string();
        let id_col_default = "<no id col>".to_string();
        let effect_col_default = "<no effect col>".to_string();
        let se_col_default = "<no se col>".to_string();
        for trait_name in &self.trait_names {
            let file = self.files.get(trait_name).unwrap_or(&file_default);
            let id_col = self.id_cols.get(trait_name).unwrap_or(&id_col_default);
            let effect_col = self.effect_cols.get(trait_name).unwrap_or(&effect_col_default);
            let se_col = self.se_cols.get(trait_name).unwrap_or(&se_col_default);
            println!("{}\t{}\t{}\t{}\t{}", trait_name, id_col, effect_col, se_col, file);
        }
    }
    fn build_mocasa_gwas_configs(&self) -> Result<Vec<GwasConfig>, Error> {
        let mut gwas_configs: Vec<GwasConfig> = Vec::new();
        let default_cols = GwasCols::default();
        for name in self.trait_names.iter() {
            let name = name.clone();
            let file = self.files.get(&name).cloned().ok_or_else(||
                Error::from(format!("No file specified for {}", name)))?;
            let id = self.id_cols.get(&name).cloned().unwrap_or(default_cols.id.clone());
            let effect =
                self.effect_cols.get(&name).cloned().unwrap_or(default_cols.effect.clone());
            let se = self.se_cols.get(&name).cloned().unwrap_or(default_cols.se.clone());
            let cols = Some(GwasCols { id, effect, se });
            gwas_configs.push(GwasConfig { name, file, cols })
        }
        Ok(gwas_configs)
    }
    fn build_mocasa_config(&self, options: &ImportPhenetOptions, phenet_opts: PhenetOpts)
                           -> Result<Config, Error> {
        let trace: Option<String> = None;
        let params = options.params_file.clone();
        let files = FilesConfig { trace, params };
        let gwas = self.build_mocasa_gwas_configs()?;
        let n_steps_burn_in = Some(defaults::shared::N_STEPS_BURN_IN);
        let n_steps_effective_burn_in = Some(defaults::shared::N_STEPS_EFFECTIVE_BURN_IN);
        let var_ratio_burn_in = Some(defaults::shared::VAR_RATIO_BURN_IN);
        let n_endos = self.endo_names.len();
        let t_pinned: Option<bool> = None;
        let shared =
            SharedConfig {
                n_endos, n_steps_burn_in, n_steps_effective_burn_in, var_ratio_burn_in, t_pinned
            };
        let PhenetOpts { var_id_file, .. } = phenet_opts;
        let ids_file = var_id_file;
        let n_samples_per_iteration = defaults::train::N_SAMPLES_PER_ITERATION;
        let n_iterations_per_round = defaults::train::N_ITERATIONS_PER_ROUND;
        let n_rounds = defaults::train::N_ROUNDS;
        let normalize_mu_to_one = true;
        let params_trace_file: Option<String> = None;
        let train =
            TrainConfig {
                ids_file,
                n_samples_per_iteration,
                n_iterations_per_round,
                n_rounds,
                normalize_mu_to_one,
                params_trace_file,
            };
        let params_override: Option<ParamsOverride> = None;
        let n_samples = Some(defaults::classify::N_SAMPLES);
        let n_samples_effective = Some(defaults::classify::N_SAMPLES_EFFECTIVE);
        let var_ratio: Option<f64> = Some(defaults::classify::VAR_RATIO);
        let n_parallel = defaults::classify::N_PARALLEL;
        let out_file = options.out_file.clone();
        let only_ids: Option<Vec<String>> = None;
        let only_ids_file: Option<String> = None;
        let trace_ids: Option<Vec<String>> = None;
        let classify =
            ClassifyConfig {
                params_override,
                n_samples,
                n_samples_effective,
                var_ratio,
                n_parallel,
                out_file,
                only_ids,
                only_ids_file,
                trace_ids,
            };
        Ok(Config { files, gwas, shared, train, classify })
    }

    fn got_some_params(&self) -> bool {
        !self.betas.is_empty() || !self.vars.is_empty() || !self.means.is_empty()
    }
    fn build_params(&self) -> Result<Params, Error> {
        let trait_names = Arc::new(self.trait_names.clone());
        let mut mus: Vec<f64> = Vec::with_capacity(self.endo_names.len());
        let mut taus: Vec<f64> = Vec::with_capacity(self.endo_names.len());
        for endo_name in &self.endo_names {
            let mu: f64 = num_or_error(&self.means, endo_name, keys::MEAN)?;
            let tau: f64 = num_or_error(&self.vars, endo_name, keys::VAR)?.sqrt();
            mus.push(mu);
            taus.push(tau);
        }
        let n_endo = self.endo_names.len();
        let n_traits = trait_names.len();
        let mut betas: Matrix = Matrix::fill(n_endo, n_traits, |_, _| 0.0);
        let mut sigmas: Vec<f64> = Vec::with_capacity(trait_names.len());
        for (i_trait, trait_name) in trait_names.iter().enumerate() {
            let beta: f64 = num_or_error(&self.betas, trait_name, keys::BETA)?;
            let sigma: f64 = num_or_error(&self.vars, trait_name, keys::VAR)?.sqrt();
            for i_endo in 0..n_endo {
                betas[i_endo][i_trait] = beta;
            }
            sigmas.push(sigma);
        }
        Ok(Params { trait_names, mus, taus, betas, sigmas })
    }
}

fn num_or_error(nums: &BTreeMap<String, f64>, name: &str, kind: &str) -> Result<f64, Error> {
    match nums.get(name) {
        None => { Err(Error::from(format!("Missing value for {} of {}.", kind, name))) }
        Some(&num) => { Ok(num) }
    }
}

pub(crate) fn import_phenet(options: &ImportPhenetOptions) -> Result<(), Error> {
    let phenet_opts = read_opts_file(&options.phenet_file)?;
    let mut config_builder = ConfigBuilder::new();
    for config_file in phenet_opts.config_files.iter() {
        config_builder.read_phenet_config(config_file)?;
    }
    if let Some(output_file) = &phenet_opts.output_file {
        config_builder.read_phenet_config_optional(output_file)?;
    }
    config_builder.report();
    let config = config_builder.build_mocasa_config(options, phenet_opts)?;
    let config_file = &options.config_file;
    write_config(&config, config_file)?;
    if config_builder.got_some_params() {
        match config_builder.build_params() {
            Ok(params) => {
                let params_string = serde_json::to_string(&params)?;
                let params_file =
                    for_file(&options.params_file, File::create(&options.params_file))?;
                let mut params_writer = BufWriter::new(params_file);
                params_writer.write_all(params_string.as_bytes())?;
            }
            Err(error) => {
                println!("Warning: no parameters written: {}", error)
            }
        }
    } else {
        println!("No parameter file written, because no params in phenet files given.")
    }
    Ok(())
}

fn missing_option_error(key: &str) -> Error {
    Error::from(format!("Missing option {}", key))
}

fn read_opts_file(opts_file: &str) -> Result<PhenetOpts, Error> {
    let mut var_id_file: Option<String> = None;
    let mut config_files: Vec<String> = Vec::new();
    let mut output_file: Option<String> = None;
    let opts_file = for_file(opts_file, File::open(opts_file))?;
    for line in BufReader::new(opts_file).lines() {
        let line = line?;
        let line = line.trim();
        if !line.starts_with('#') && !line.is_empty() {
            let mut parts_iter = line.split(char::is_whitespace);
            match (parts_iter.next(), parts_iter.next(), parts_iter.next()) {
                (Some(key), Some(value), None) => {
                    match key {
                        keys::VAR_ID_FILE => { var_id_file = Some(value.to_string()) }
                        keys::CONFIG_FILE => { config_files.push(value.to_string()) }
                        keys::OUTPUT_FILE => { output_file = Some(value.to_string()) }
                        _ => { println!("Ignoring option: {}: {}", key, value) }
                    }
                }
                _ => {
                    Err(Error::from(
                        format!("Expected key/value pair, but got '{}'.", line)
                    ))?;
                }
            }
        }
    }
    let var_id_file =
        match var_id_file {
            None => {
                println!("Warning: no variant id file ({}) specified.", keys::VAR_ID_FILE);
                "".to_string()
            }
            Some(var_id_file) => { var_id_file }
        };
    // var_id_file.ok_or_else(|| missing_option_error(keys::VAR_ID_FILE))?;
    if config_files.is_empty() {
        Err(missing_option_error(keys::CONFIG_FILE))?;
    }
    Ok(PhenetOpts { var_id_file, config_files, output_file })
}