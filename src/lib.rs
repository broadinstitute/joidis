use crate::error::Error;
use crate::options::cli::get_cli_options;
use crate::options::config::load_config;

mod model;
mod observables;
mod options;
mod error;

pub fn run() -> Result<(), Error> {
    let cli_options = get_cli_options()?;
    let config = load_config(&cli_options.config_file)?;
    println!("Number of GWAS files: {}", config.gwas.len());
    Ok(())
}
