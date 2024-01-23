use crate::train::params::Params;

pub(crate) fn calculate_mu(params: &Params, betas: &[f64], ses: &[f64]) -> f64 {
    let numerator: f64 =
        params.betas.iter().zip(params.sigmas.iter()).zip(betas.iter())
            .zip(ses.iter())
            .map(|(((&beta, &sigma), &o), &se)| {
                beta * o / (sigma.powi(2) + se.powi(2))
            }).sum();
    let denominator: f64 =
        params.betas.iter().zip(params.sigmas.iter()).zip(ses.iter())
            .map(|((&beta, &sigma), &se)| {
                beta.powi(2) / (sigma.powi(2) + se.powi(2))
            }).sum();
    numerator / denominator
}