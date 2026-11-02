// SPDX-License-Identifier: GPL-3.0-or-later
//! Experimental observed-pair likelihood for supplied response priors.
//!
//! This evaluates a research objective, not a qualified hyperparameter learner.
//! Selection conditioning uses the JOINT response prior, not the null prior.
use crate::observation::{PairTable, SampleConfiguration};
use crate::response::{JointPrior, ResponseError, ResponseModel};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseObjective {
    pub log_observation_probability: f64,
    pub unconditioned_negative_log_likelihood: f64,
    pub conditional_negative_log_likelihood: f64,
    pub observations: u64,
}

fn log_sum(mut values: Vec<f64>) -> Result<f64, ResponseError> {
    if values.iter().any(|x| x.is_nan() || *x == f64::INFINITY) {
        return Err(ResponseError::NumericalRange);
    }
    let maximum = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if maximum == f64::NEG_INFINITY {
        return Err(ResponseError::AllImpossible);
    }
    for x in &mut values {
        *x = (*x - maximum).exp();
    }
    values.sort_unstable_by(f64::total_cmp);
    Ok(maximum + values.into_iter().sum::<f64>().ln())
}

/// Evaluate -sum multiplicity*log(Q(n1,n2)) + Nobs*log(1-Q(0,0)).
/// The ordinary posterior precision checks apply; no automatic optimizer,
/// parameter uncertainty, model selection or biological calibration is implied.
pub fn evaluate(
    prior: &JointPrior,
    samples: [SampleConfiguration; 2],
    table: &PairTable,
) -> Result<ResponseObjective, ResponseError> {
    if table.observations() > 1_u64 << 53 {
        return Err(ResponseError::NumericalRange);
    }
    let model = ResponseModel::new(prior.clone(), samples[0], samples[1], 0.0, 1.0)?;
    let mut observed = Vec::new();
    observed
        .try_reserve_exact(prior.states().len())
        .map_err(|_| ResponseError::AllocationFailure)?;
    for state in prior.states() {
        let zero = samples[0]
            .acquisition
            .log_pmf(0, state.frequency1, samples[0].depth)?
            + samples[1]
                .acquisition
                .log_pmf(0, state.frequency2, samples[1].depth)?;
        // Sum P(observed|state) directly to retain accuracy when Q00 is near1.
        let log_observed = (-zero.exp_m1()).ln();
        observed.push(state.mass.ln() - prior.raw_mass().ln() + log_observed);
    }
    let log_observed = log_sum(observed)?;
    if log_observed > 64.0 * f64::EPSILON {
        return Err(ResponseError::InvalidPosterior);
    }
    let log_observed = log_observed.min(0.0);
    let mut nll = 0.0;
    for row in table.records() {
        nll -= row.multiplicity as f64 * model.infer(row.count1, row.count2)?.log_evidence;
    }
    let conditional = nll + table.observations() as f64 * log_observed;
    if !nll.is_finite() || !conditional.is_finite() {
        return Err(ResponseError::NumericalRange);
    }
    Ok(ResponseObjective {
        log_observation_probability: log_observed,
        unconditioned_negative_log_likelihood: nll,
        conditional_negative_log_likelihood: conditional,
        observations: table.observations(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::acquisition::Acquisition;
    use crate::observation::ObservedCountPair;
    use crate::response::JointState;

    #[test]
    fn joint_selection_normalizes_observed_distribution_and_objective() {
        let prior = JointPrior::new(
            &[
                JointState {
                    frequency1: 0.1,
                    frequency2: 0.1,
                    mass: 0.6,
                },
                JointState {
                    frequency1: 0.1,
                    frequency2: 0.5,
                    mass: 0.2,
                },
                JointState {
                    frequency1: 0.5,
                    frequency2: 0.1,
                    mass: 0.2,
                },
            ],
            true,
        )
        .unwrap();
        let samples = [SampleConfiguration {
            acquisition: Acquisition::Poisson,
            depth: 2.0,
        }; 2];
        let table = PairTable::new(&[
            ObservedCountPair {
                count1: 1,
                count2: 0,
                multiplicity: 3,
            },
            ObservedCountPair {
                count1: 0,
                count2: 2,
                multiplicity: 2,
            },
        ])
        .unwrap();
        let objective = evaluate(&prior, samples, &table).unwrap();
        let denominator = 1.0 - (0.6 * (-0.4_f64).exp() + 0.4 * (-1.2_f64).exp());
        assert!((objective.log_observation_probability - denominator.ln()).abs() < 1e-15);
        assert!(
            (objective.conditional_negative_log_likelihood
                - objective.unconditioned_negative_log_likelihood
                - 5.0 * denominator.ln())
            .abs()
                < 2e-15
        );
        let model = ResponseModel::new(prior, samples[0], samples[1], 0.0, 1.0).unwrap();
        let mut total = 0.0;
        for a in 0..24 {
            for b in 0..24 {
                if a + b > 0 {
                    total += (model.infer(a, b).unwrap().log_evidence - denominator.ln()).exp();
                }
            }
        }
        assert!((total - 1.0).abs() < 3e-15);
    }

    #[test]
    fn rare_observation_and_raw_mass_roundoff_remain_finite() {
        let prior = JointPrior::new(
            &[JointState {
                frequency1: 0.1,
                frequency2: 0.2,
                mass: 1.0 + 2e-13,
            }],
            false,
        )
        .unwrap();
        let table = PairTable::new(&[ObservedCountPair {
            count1: 1,
            count2: 0,
            multiplicity: 1,
        }])
        .unwrap();
        for depth in [1e-100, 1000.0] {
            let samples = [SampleConfiguration {
                acquisition: Acquisition::Poisson,
                depth,
            }; 2];
            let result = evaluate(&prior, samples, &table).unwrap();
            let expected = (-(-0.3 * depth).exp_m1()).ln();
            assert!((result.log_observation_probability - expected).abs() < 1e-12);
            assert!(result.conditional_negative_log_likelihood.is_finite());
        }
    }
}
