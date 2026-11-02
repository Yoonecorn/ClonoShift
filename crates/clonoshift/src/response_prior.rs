// SPDX-License-Identifier: GPL-3.0-or-later
//! Scout midpoint prior. No fitted-noise coupling or hyperparameter estimation.
use crate::numerics::FrequencyPrior;
use crate::observation::Grid;
use crate::response::{JointPrior, JointState, ResponseError};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Baseline {
    pub exponent: f64,
    pub minimum: f64,
}
#[derive(Debug, Clone, Copy)]
pub struct ResponsePriorConfig {
    pub response_probability: f64,
    pub scale: f64,
    pub points: usize,
    pub stride: usize,
    pub maximum_shift: usize,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PriorBuildError {
    InvalidConfiguration,
    IncompatibleBaselines,
    NumericalRange,
    AllocationFailure,
    Posterior(ResponseError),
}
#[derive(Debug, Clone)]
pub struct Slice {
    pub shift: usize,
    pub response: f64,
    /// One signed slice (the zero slice is counted once).
    pub marginal_mass: f64,
    pub log_midpoint_normalizer: f64,
    pub midpoint_nodes: usize,
    pub singleton: bool,
}
#[derive(Debug, Clone)]
pub struct BuiltPrior {
    pub prior: JointPrior,
    pub grid: Grid,
    pub slices: Vec<Slice>,
}
fn reserve<T>(n: usize) -> Result<Vec<T>, PriorBuildError> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| PriorBuildError::AllocationFailure)?;
    Ok(v)
}
fn lse(v: &[f64]) -> f64 {
    let m = v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    m + v.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
}
/// Dual baselines are an explicit compatibility check, never a synthesis rule.
/// Discrete slab node masses are proportional to exp(-abs(r)/scale), including zero.
/// Validate geometry even for a diagonal-only prior. Singleton slices use point mass one.
pub fn build_symmetric_prior(
    baseline: Baseline,
    second_baseline: Option<Baseline>,
    config: ResponsePriorConfig,
) -> Result<BuiltPrior, PriorBuildError> {
    use PriorBuildError::*;
    if second_baseline.is_some_and(|other| other != baseline) {
        return Err(IncompatibleBaselines);
    }
    let c = config;
    let last = c
        .stride
        .checked_mul(c.maximum_shift)
        .ok_or(InvalidConfiguration)?;
    if !c.response_probability.is_finite()
        || !(0.0..=1.0).contains(&c.response_probability)
        || !c.scale.is_finite()
        || c.scale <= 0.0
        || c.stride == 0
        || c.points < 2
        || c.points > 1025
        || last.checked_mul(2).ok_or(InvalidConfiguration)? >= c.points
        || (c.response_probability > 0.0 && c.maximum_shift == 0)
    {
        return Err(InvalidConfiguration);
    }
    // Bound the unchanged posterior constructor's quadratic validation cost.
    let state_count = if c.response_probability == 0.0 {
        c.points
    } else {
        (0..=c.maximum_shift)
            .map(|j| (c.points - 2 * j * c.stride) * if j == 0 { 1 } else { 2 })
            .sum()
    };
    if state_count > 40_000 {
        return Err(InvalidConfiguration);
    }
    let density = FrequencyPrior::new(baseline.exponent, baseline.minimum)
        .map_err(|_| InvalidConfiguration)?;
    let grid = Grid::log_uniform(baseline.minimum, c.points).map_err(|_| InvalidConfiguration)?;
    let f = grid.frequencies();
    let dx = -baseline.minimum.ln() / (c.points - 1) as f64;
    let mut slab = reserve(c.maximum_shift + 1)?;
    for k in 0..=c.maximum_shift {
        slab.push(-2.0 * (k * c.stride) as f64 * dx / c.scale);
    }
    // The zero node occurs once; every positive node occurs twice.
    let mut signed = reserve(2 * c.maximum_shift + 1)?;
    signed.push(slab[0]);
    for &v in &slab[1..] {
        signed.extend([v, v]);
    }
    let norm = lse(&signed);
    let mut states = reserve(c.points * (2 * c.maximum_shift + 1))?;
    let mut slices = reserve(c.maximum_shift + 1)?;
    for (j, &slab_log) in slab.iter().enumerate() {
        let k = j * c.stride;
        let marginal = if j == 0 {
            1.0 - c.response_probability + c.response_probability * (slab_log - norm).exp()
        } else {
            c.response_probability * (slab_log - norm).exp()
        };
        if marginal == 0.0 {
            if c.response_probability > 0.0 {
                return Err(NumericalRange);
            }
            continue;
        }
        let nodes = c.points - 2 * k;
        let mut logs = reserve(nodes)?;
        for i in k..c.points - k {
            let weight = if nodes == 1 {
                1.0
            } else if i == k {
                (f[i + 1].ln() - f[i].ln()) / 2.0
            } else if i == c.points - k - 1 {
                (f[i].ln() - f[i - 1].ln()) / 2.0
            } else {
                (f[i + 1].ln() - f[i - 1].ln()) / 2.0
            };
            logs.push(
                density.log_density(f[i]).map_err(|_| NumericalRange)? + f[i].ln() + weight.ln(),
            );
        }
        let z = lse(&logs);
        if !z.is_finite() {
            return Err(NumericalRange);
        }
        for (offset, &log) in logs.iter().enumerate() {
            let i = k + offset;
            let mass = (marginal.ln() + log - z).exp();
            if !mass.is_finite() || mass <= 0.0 {
                return Err(NumericalRange);
            }
            let state = JointState {
                frequency1: f[i - k],
                frequency2: f[i + k],
                mass,
            };
            states.push(state);
            if k > 0 {
                states.push(JointState {
                    frequency1: state.frequency2,
                    frequency2: state.frequency1,
                    mass,
                });
            }
        }
        slices.push(Slice {
            shift: k,
            response: 2.0 * k as f64 * dx,
            marginal_mass: marginal,
            log_midpoint_normalizer: z,
            midpoint_nodes: nodes,
            singleton: nodes == 1,
        });
    }
    let prior = JointPrior::new(&states, true).map_err(Posterior)?;
    Ok(BuiltPrior {
        prior,
        grid,
        slices,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn impossible_allocation_is_a_classified_error() {
        assert_eq!(
            reserve::<JointState>(usize::MAX).unwrap_err(),
            PriorBuildError::AllocationFailure
        );
    }
}
