use crate::data::gwas::GwasCols;
use crate::error::{for_file, Error};
use crate::options::cli::ImportPhenetOptions;
use crate::options::config::{ClassifyConfig, Config, FilesConfig, GwasConfig, TrainConfig};
use crate::params::{Params, ParamsOverride};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;

mod defaults {
    pub(crate) mod train {
        pub(crate) const N_STEPS_BURN_IN: usize = 10000;
        pub(crate) const N_SAMPLES_PER_ITERATION: usize = 100;
        pub(crate) const N_ITERATIONS_PER_ROUND: usize = 1000;
        pub(crate) const N_ROUNDS: usize = 10000;
    }

    pub(crate) mod classify {
        pub(crate) const N_STEPS_BURN_IN: usize = 10000;
        pub(crate) const N_SAMPLES: usize = 100000;
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
    endo_name: Option<String>,
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
        let endo_name: Option<String> = None;
        let files: BTreeMap<String, String> = BTreeMap::new();
        let id_cols: BTreeMap<String, String> = BTreeMap::new();
        let effect_cols: BTreeMap<String, String> = BTreeMap::new();
        let se_cols: BTreeMap<String, String> = BTreeMap::new();
        let betas: BTreeMap<String, f64> = BTreeMap::new();
        let vars: BTreeMap<String, f64> = BTreeMap::new();
        let means: BTreeMap<String, f64> = BTreeMap::new();
        ConfigBuilder {
            trait_names,
            endo_name,
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
                match (
                    parts_iter.next(),
                    parts_iter.next(),
                    parts_iter.next(),
                    parts_iter.next(),
                ) {
                    (Some(part1), Some(part2), Some(part3), None) => match (part1, part2, part3) {
                        (trait_name, keys::DECLARE, keys::TRAIT) => {
                            let trait_name = trait_name.to_string();
                            if !self.trait_names.contains(&trait_name) {
                                self.trait_names.push(trait_name.to_string());
                            }
                        }
                        (endo_name, keys::DECLARE, keys::ENDO) => {
                            self.endo_name = Some(endo_name.to_string());
                        }
                        (trait_name, keys::FILE, file) => {
                            self.files.insert(trait_name.to_string(), file.to_string());
                        }
                        (trait_name, keys::ID_COL, id_col) => {
                            self.id_cols
                                .insert(trait_name.to_string(), id_col.to_string());
                        }
                        (trait_name, keys::EFFECT_COL, effect_col) => {
                            self.effect_cols
                                .insert(trait_name.to_string(), effect_col.to_string());
                        }
                        (trait_name, keys::SE_COL, se_col) => {
                            self.se_cols
                                .insert(trait_name.to_string(), se_col.to_string());
                        }
                        (trait_name, keys::BETA, beta) => {
                            self.betas
                                .insert(trait_name.to_string(), beta.parse::<f64>()?);
                        }
                        (trait_name, keys::VAR, var) => {
                            self.vars
                                .insert(trait_name.to_string(), var.parse::<f64>()?);
                        }
                        (trait_name, keys::MEAN, mean) => {
                            self.means
                                .insert(trait_name.to_string(), mean.parse::<f64>()?);
                        }
                        _ => {
                            println!("Ignoring line '{}'.", line)
                        }
                    },
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
        if let Some(endo_name) = &self.endo_name {
            println!("Endo name: {}", endo_name)
        }
        let file_default = "<no file>".to_string();
        let id_col_default = "<no id col>".to_string();
        let effect_col_default = "<no effect col>".to_string();
        let se_col_default = "<no se col>".to_string();
        for trait_name in &self.trait_names {
            let file = self.files.get(trait_name).unwrap_or(&file_default);
            let id_col = self.id_cols.get(trait_name).unwrap_or(&id_col_default);
            let effect_col = self
                .effect_cols
                .get(trait_name)
                .unwrap_or(&effect_col_default);
            let se_col = self.se_cols.get(trait_name).unwrap_or(&se_col_default);
            println!(
                "{}\t{}\t{}\t{}\t{}",
                trait_name, id_col, effect_col, se_col, file
            );
        }
    }
    fn build_mocasa_gwas_configs(&self) -> Result<Vec<GwasConfig>, Error> {
        let mut gwas_configs: Vec<GwasConfig> = Vec::new();
        let default_cols = GwasCols::default();
        for name in self.trait_names.iter() {
            let name = name.clone();
            let file = self
                .files
                .get(&name)
                .cloned()
                .ok_or_else(|| Error::from(format!("No file specified for {}", name)))?;
            let id = self
                .id_cols
                .get(&name)
                .cloned()
                .unwrap_or(default_cols.id.clone());
            let effect = self
                .effect_cols
                .get(&name)
                .cloned()
                .unwrap_or(default_cols.effect.clone());
            let se = self
                .se_cols
                .get(&name)
                .cloned()
                .unwrap_or(default_cols.se.clone());
            let cols = Some(GwasCols { id, effect, se });
            gwas_configs.push(GwasConfig { name, file, cols })
        }
        Ok(gwas_configs)
    }
    fn build_mocasa_config(
        &self,
        options: &ImportPhenetOptions,
        phenet_opts: PhenetOpts,
    ) -> Result<Config, Error> {
        let trace: Option<String> = None;
        let params = options.params_file.clone();
        let files = FilesConfig { trace, params };
        let gwas = self.build_mocasa_gwas_configs()?;
        let PhenetOpts { var_id_file, .. } = phenet_opts;
        let ids_file = var_id_file;
        let n_steps_burn_in = defaults::train::N_STEPS_BURN_IN;
        let n_samples_per_iteration = defaults::train::N_SAMPLES_PER_ITERATION;
        let n_iterations_per_round = defaults::train::N_ITERATIONS_PER_ROUND;
        let n_rounds = defaults::train::N_ROUNDS;
        let normalize_mu_to_one = true;
        let params_trace_file: Option<String> = None;
        let t_pinned: Option<bool> = None;
        let train = TrainConfig {
            ids_file,
            n_steps_burn_in,
            n_samples_per_iteration,
            n_iterations_per_round,
            n_rounds,
            normalize_mu_to_one,
            params_trace_file,
            t_pinned,
        };
        let params_override: Option<ParamsOverride> = None;
        let n_steps_burn_in = defaults::classify::N_STEPS_BURN_IN;
        let n_samples = defaults::classify::N_SAMPLES;
        let out_file = options.out_file.clone();
        let trace_ids: Option<Vec<String>> = None;
        let classify = ClassifyConfig {
            params_override,
            n_steps_burn_in,
            n_samples,
            out_file,
            trace_ids,
            t_pinned,
        };
        Ok(Config {
            files,
            gwas,
            train,
            classify,
        })
    }

    fn got_some_params(&self) -> bool {
        !self.betas.is_empty() || !self.vars.is_empty() || !self.means.is_empty()
    }
    fn build_params(&self) -> Result<Params, Error> {
        let endo = self
            .endo_name
            .as_ref()
            .ok_or_else(|| Error::from("No name for endo phenotype specified."))?;
        let trait_names = Arc::new(self.trait_names.clone());
        let mu: f64 = num_or_error(&self.means, endo, keys::MEAN)?;
        let tau: f64 = num_or_error(&self.vars, endo, keys::VAR)?.sqrt();
        let mut betas: Vec<f64> = Vec::with_capacity(trait_names.len());
        let mut sigmas: Vec<f64> = Vec::with_capacity(trait_names.len());
        for trait_name in trait_names.iter() {
            let beta: f64 = num_or_error(&self.betas, trait_name, keys::BETA)?;
            let sigma: f64 = num_or_error(&self.vars, trait_name, keys::VAR)?.sqrt();
            betas.push(beta);
            sigmas.push(sigma);
        }
        Ok(Params {
            trait_names,
            mu,
            tau,
            betas,
            sigmas,
        })
    }
}

fn num_or_error(nums: &BTreeMap<String, f64>, name: &str, kind: &str) -> Result<f64, Error> {
    match nums.get(name) {
        None => Err(Error::from(format!(
            "Missing value for {} of {}.",
            kind, name
        ))),
        Some(&num) => Ok(num),
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
    let config_string = toml::to_string(&config)?;
    let config_file = for_file(&options.config_file, File::create(&options.config_file))?;
    let mut writer = BufWriter::new(config_file);
    writer.write_all(config_string.as_bytes())?;
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
                (Some(key), Some(value), None) => match key {
                    keys::VAR_ID_FILE => var_id_file = Some(value.to_string()),
                    keys::CONFIG_FILE => config_files.push(value.to_string()),
                    keys::OUTPUT_FILE => output_file = Some(value.to_string()),
                    _ => {
                        println!("Ignoring option: {}: {}", key, value)
                    }
                },
                _ => {
                    Err(Error::from(format!(
                        "Expected key/value pair, but got '{}'.",
                        line
                    )))?;
                }
            }
        }
    }
    let var_id_file = match var_id_file {
        None => {
            println!(
                "Warning: no variant id file ({}) specified.",
                keys::VAR_ID_FILE
            );
            "".to_string()
        }
        Some(var_id_file) => var_id_file,
    };
    // var_id_file.ok_or_else(|| missing_option_error(keys::VAR_ID_FILE))?;
    if config_files.is_empty() {
        Err(missing_option_error(keys::CONFIG_FILE))?;
    }
    Ok(PhenetOpts {
        var_id_file,
        config_files,
        output_file,
    })
}
