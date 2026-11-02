// SPDX-License-Identifier: GPL-3.0-or-later
//! Explicit two-sample observation conditioning and fixed-parameter objectives.
//! Scientific and numerical contracts: docs/METHODS.md.
use crate::acquisition::{Acquisition, AcquisitionError};
use crate::numerics::{FrequencyPrior, NumericalError};
use std::fmt;

const ROUND: f64 = 64.0 * f64::EPSILON;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationError {
    InvalidGrid,
    DimensionMismatch,
    AllocationFailure,
    InvalidObservation,
    InvalidMultiplicity,
    CountOverflow,
    EmptyObservations,
    InvalidP00,
    InvalidDenominator,
    UnexpectedZeroProbability,
    InvalidProbability,
    InconsistentNormalization,
    NumericalRange,
    Acquisition(AcquisitionError),
    Prior(NumericalError),
}
impl fmt::Display for ObservationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "observation error: {self:?}")
    }
}
impl std::error::Error for ObservationError {}
impl From<AcquisitionError> for ObservationError {
    fn from(e: AcquisitionError) -> Self {
        Self::Acquisition(e)
    }
}
impl From<NumericalError> for ObservationError {
    fn from(e: NumericalError) -> Self {
        Self::Prior(e)
    }
}
type Result<T> = std::result::Result<T, ObservationError>;
fn storage<T>(n: usize) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| ObservationError::AllocationFailure)?;
    Ok(values)
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(ObservationError::NumericalRange)
    }
}
fn lse<I: Iterator<Item = f64> + Clone>(values: I) -> Result<f64> {
    let mut max = f64::NEG_INFINITY;
    for value in values.clone() {
        if value.is_nan() || value == f64::INFINITY {
            return Err(ObservationError::NumericalRange);
        }
        max = max.max(value);
    }
    if max == f64::NEG_INFINITY {
        return Ok(max);
    }
    finite(max + values.map(|v| (v - max).exp()).sum::<f64>().ln())
}
fn log_add(a: f64, b: f64) -> Result<f64> {
    lse([a, b].into_iter())
}

// Return log integral and log posterior mean f. Center acquisition logs before
// adding prior weights, then center again before adding log f. Subtracting two
// enormous final log integrals would erase the O(1) posterior moment ratio.
fn pair_integrals(raw: &[f64], logs: &[f64], log_fs: &[f64], log_mass: f64) -> Result<(f64, f64)> {
    if raw.len() != logs.len() || raw.len() != log_fs.len() {
        return Err(ObservationError::DimensionMismatch);
    }
    let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max == f64::NEG_INFINITY {
        return Ok((max, max));
    }
    let mut centered = storage(raw.len())?;
    for (&r, &p) in raw.iter().zip(logs) {
        centered.push(if p == f64::NEG_INFINITY {
            p
        } else {
            finite(r + (p - max))?
        });
    }
    let offset = centered.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    for value in &mut centered {
        *value -= offset;
    }
    let normalizer = lse(centered.iter().copied())?;
    let mean = finite(lse(centered.iter().zip(log_fs).map(|(r, f)| r + f))? - normalizer)?;
    let q = finite(max + finite(offset + normalizer - log_mass)?)?;
    Ok((q, mean))
}

/// Trapezoids in x=ln(f), with strictly increasing explicit frequencies.
#[derive(Debug, Clone)]
pub struct Grid {
    frequencies: Vec<f64>,
    log_frequencies: Vec<f64>,
    weights: Vec<f64>,
}
impl Grid {
    pub fn from_frequencies(frequencies: &[f64]) -> Result<Self> {
        if frequencies.len() < 2 || frequencies.last() != Some(&1.0) {
            return Err(ObservationError::InvalidGrid);
        }
        let mut fs = storage(frequencies.len())?;
        let mut xs = storage(frequencies.len())?;
        let mut last = 0.0;
        for &f in frequencies {
            if !f.is_finite() || f <= last || f > 1.0 {
                return Err(ObservationError::InvalidGrid);
            }
            fs.push(f);
            xs.push(f.ln());
            last = f;
        }
        let mut weights = storage(fs.len())?;
        for i in 0..fs.len() {
            let left = if i == 0 { xs[i] } else { xs[i - 1] };
            let right = if i + 1 == fs.len() { xs[i] } else { xs[i + 1] };
            let w = (right - left) / 2.0;
            if !w.is_finite() || w <= 0.0 {
                return Err(ObservationError::InvalidGrid);
            }
            weights.push(w);
        }
        Ok(Self {
            frequencies: fs,
            log_frequencies: xs,
            weights,
        })
    }
    /// Supplied weights must exactly match this platform's documented f64 rule.
    /// Neither category/node ordering nor caller-provided values are repaired.
    pub fn with_trapezoid_weights(frequencies: &[f64], weights: &[f64]) -> Result<Self> {
        if frequencies.len() != weights.len() {
            return Err(ObservationError::DimensionMismatch);
        }
        let grid = Self::from_frequencies(frequencies)?;
        if grid.weights != weights {
            return Err(ObservationError::InvalidGrid);
        }
        Ok(grid)
    }
    pub fn log_uniform(minimum: f64, points: usize) -> Result<Self> {
        if !minimum.is_finite() || minimum <= 0.0 || minimum >= 1.0 || points < 2 {
            return Err(ObservationError::InvalidGrid);
        }
        let mut fs = storage(points)?;
        for i in 0..points {
            fs.push(if i == 0 {
                minimum
            } else if i + 1 == points {
                1.0
            } else {
                (minimum.ln() * (1.0 - i as f64 / (points - 1) as f64)).exp()
            });
        }
        Self::from_frequencies(&fs)
    }
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
}
#[derive(Debug, Clone, Copy)]
pub struct SampleConfiguration {
    pub acquisition: Acquisition,
    pub depth: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LogProbability {
    log_value: f64,
    roundoff_corrections: u64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirectProbability {
    pub value: f64,
    pub underflowed: bool,
}
impl LogProbability {
    pub fn log_value(self) -> f64 {
        self.log_value
    }
    pub fn roundoff_corrections(self) -> u64 {
        self.roundoff_corrections
    }
    pub fn direct(self) -> DirectProbability {
        let value = self.log_value.exp();
        DirectProbability {
            value,
            underflowed: value == 0.0 && self.log_value.is_finite(),
        }
    }
}
fn probability(value: f64, zero_allowed: bool) -> Result<LogProbability> {
    if value == f64::NEG_INFINITY {
        return if zero_allowed {
            Ok(LogProbability {
                log_value: value,
                roundoff_corrections: 0,
            })
        } else {
            Err(ObservationError::UnexpectedZeroProbability)
        };
    }
    if !value.is_finite() || value > ROUND {
        return Err(ObservationError::InvalidProbability);
    }
    Ok(LogProbability {
        log_value: if value > 0.0 { 0.0 } else { value },
        roundoff_corrections: u64::from(value > 0.0),
    })
}
fn valid_denominator(log_s: f64) -> Result<f64> {
    if !log_s.is_finite() || log_s > 0.0 {
        Err(ObservationError::InvalidDenominator)
    } else {
        Ok(log_s)
    }
}
fn valid_mass(log_mass: f64) -> Result<()> {
    if !log_mass.is_finite() || log_mass.abs() > ROUND {
        Err(ObservationError::InconsistentNormalization)
    } else {
        Ok(())
    }
}
fn conditioning(log_p00: f64) -> Result<(f64, f64)> {
    if !log_p00.is_finite() || log_p00 >= 0.0 {
        return Err(ObservationError::InvalidP00);
    }
    let p00 = log_p00.exp();
    if !(0.0..1.0).contains(&p00) {
        return Err(ObservationError::InvalidP00);
    }
    let log_s = (-log_p00.exp_m1()).ln();
    Ok((p00, valid_denominator(log_s)?))
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObservedCountPair {
    pub count1: u64,
    pub count2: u64,
    pub multiplicity: u64,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairTable {
    records: Vec<ObservedCountPair>,
    observations: u64,
}
impl PairTable {
    pub fn new(records: &[ObservedCountPair]) -> Result<Self> {
        let mut copy = storage(records.len())?;
        let mut observations = 0_u64;
        for &r in records {
            if r.count1 == 0 && r.count2 == 0 {
                return Err(ObservationError::InvalidObservation);
            }
            if r.multiplicity == 0 {
                return Err(ObservationError::InvalidMultiplicity);
            }
            observations = observations
                .checked_add(r.multiplicity)
                .ok_or(ObservationError::CountOverflow)?;
            copy.push(r);
        }
        copy.sort_unstable_by_key(|r| (r.count1, r.count2));
        let mut grouped: Vec<ObservedCountPair> = storage(copy.len())?;
        for r in copy {
            if let Some(previous) = grouped.last_mut() {
                if (previous.count1, previous.count2) == (r.count1, r.count2) {
                    previous.multiplicity = previous
                        .multiplicity
                        .checked_add(r.multiplicity)
                        .ok_or(ObservationError::CountOverflow)?;
                    continue;
                }
            }
            grouped.push(r);
        }
        Ok(Self {
            records: grouped,
            observations,
        })
    }
    pub fn records(&self) -> &[ObservedCountPair] {
        &self.records
    }
    pub fn observations(&self) -> u64 {
        self.observations
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairEvaluation {
    pub unconditioned: LogProbability,
    pub conditioned: LogProbability,
    pub log_posterior_mean_frequency: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObjectiveSummary {
    pub negative_log_likelihood: f64,
    pub mean_negative_log_likelihood: f64,
    pub c1: f64,
    pub c2: f64,
    pub log_population_count: f64,
    pub log_p00: f64,
    pub p00: f64,
    pub log_observation_probability: f64,
    pub raw_prior_log_mass: f64,
    pub observations: u64,
    pub distinct_pairs: usize,
    pub grid_points: usize,
    pub roundoff_corrections: u64,
    pub direct_probability_underflows: u64,
}
/// Fixed parameters only: no optimizer state, penalty terms or model selection.
#[derive(Debug, Clone)]
pub struct ConditionedModel {
    grid: Grid,
    sample1: SampleConfiguration,
    sample2: SampleConfiguration,
    raw_log_masses: Vec<f64>,
    log_mass: f64,
    log_p00: f64,
    p00: f64,
    log_s: f64,
    log_i00: f64,
    log_mean_frequency: f64,
}
impl ConditionedModel {
    pub fn new(
        grid: Grid,
        exponent: f64,
        sample1: SampleConfiguration,
        sample2: SampleConfiguration,
    ) -> Result<Self> {
        let prior = FrequencyPrior::new(exponent, grid.frequencies[0])?;
        let mut raw = storage(grid.frequencies.len())?;
        let mut zero_logs = storage(grid.frequencies.len())?;
        for (i, &f) in grid.frequencies.iter().enumerate() {
            raw.push(finite(
                grid.weights[i].ln() + prior.log_density(f)? + grid.log_frequencies[i],
            )?);
            let a = sample1.acquisition.log_pmf(0, f, sample1.depth)?;
            let b = sample2.acquisition.log_pmf(0, f, sample2.depth)?;
            zero_logs.push(finite(a + b)?);
        }
        let log_mass = lse(raw.iter().copied())?;
        let normalized_mass = lse(raw.iter().map(|r| r - log_mass))?;
        valid_mass(normalized_mass)?;
        let (log_p00, log_mean00) =
            pair_integrals(&raw, &zero_logs, &grid.log_frequencies, log_mass)?;
        let (p00, log_s) = conditioning(log_p00)?;
        let log_i00 = finite(log_p00 + log_mean00)?;
        let log_mean_frequency =
            finite(lse(raw.iter().zip(&grid.log_frequencies).map(|(r, f)| r + f))? - log_mass)?;
        Ok(Self {
            grid,
            sample1,
            sample2,
            raw_log_masses: raw,
            log_mass,
            log_p00,
            p00,
            log_s,
            log_i00,
            log_mean_frequency,
        })
    }
    pub fn grid(&self) -> &Grid {
        &self.grid
    }
    pub fn log_p00(&self) -> f64 {
        self.log_p00
    }
    pub fn p00(&self) -> f64 {
        self.p00
    }
    pub fn log_observation_probability(&self) -> f64 {
        self.log_s
    }
    pub fn raw_prior_log_mass(&self) -> f64 {
        self.log_mass
    }
    pub fn normalized_log_masses(&self) -> Result<Vec<f64>> {
        let mut masses = storage(self.raw_log_masses.len())?;
        masses.extend(self.raw_log_masses.iter().map(|x| x - self.log_mass));
        Ok(masses)
    }
    fn pair_moments(&self, k1: u64, k2: u64) -> Result<(LogProbability, f64)> {
        let mut logs = storage(self.grid.frequencies.len())?;
        for &f in &self.grid.frequencies {
            let a = self
                .sample1
                .acquisition
                .log_pmf(k1, f, self.sample1.depth)?;
            let b = self
                .sample2
                .acquisition
                .log_pmf(k2, f, self.sample2.depth)?;
            if a == f64::NEG_INFINITY || b == f64::NEG_INFINITY {
                logs.push(f64::NEG_INFINITY);
            } else {
                logs.push(finite(a + b)?);
            }
        }
        let (log_q, log_mean) = pair_integrals(
            &self.raw_log_masses,
            &logs,
            &self.grid.log_frequencies,
            self.log_mass,
        )?;
        Ok((probability(log_q, true)?, log_mean))
    }

    /// Includes explicit00; mathematically impossible pairs have log -infinity.
    pub fn log_pair(&self, k1: u64, k2: u64) -> Result<LogProbability> {
        Ok(self.pair_moments(k1, k2)?.0)
    }
    pub fn observed_pair(&self, k1: u64, k2: u64) -> Result<PairEvaluation> {
        if k1 == 0 && k2 == 0 {
            return Err(ObservationError::InvalidObservation);
        }
        let (unconditioned, mean) = self.pair_moments(k1, k2)?;
        if unconditioned.log_value == f64::NEG_INFINITY {
            return Err(ObservationError::UnexpectedZeroProbability);
        }
        let conditioned = probability(finite(unconditioned.log_value - self.log_s)?, false)?;
        let mean = finite(mean)?;
        if mean < self.grid.log_frequencies[0] - ROUND || mean > ROUND {
            return Err(ObservationError::InconsistentNormalization);
        }
        Ok(PairEvaluation {
            unconditioned,
            conditioned,
            log_posterior_mean_frequency: mean,
        })
    }
    pub fn evaluate(&self, table: &PairTable) -> Result<ObjectiveSummary> {
        if table.observations == 0 {
            return Err(ObservationError::EmptyObservations);
        }
        let log_population_count = (table.observations as f64).ln() - self.log_s;
        let mut nll = 0.0;
        let mut observed_moments = storage(table.records.len())?;
        let mut corrections = 0_u64;
        let mut underflows = u64::from(self.p00 == 0.0 && self.log_p00.is_finite());
        for r in &table.records {
            let pair = self.observed_pair(r.count1, r.count2)?;
            nll = finite(nll - (r.multiplicity as f64) * pair.conditioned.log_value)?;
            observed_moments.push((r.multiplicity as f64).ln() + pair.log_posterior_mean_frequency);
            corrections = corrections
                .checked_add(pair.unconditioned.roundoff_corrections)
                .and_then(|v| v.checked_add(pair.conditioned.roundoff_corrections))
                .ok_or(ObservationError::CountOverflow)?;
            underflows = underflows
                .checked_add(u64::from(pair.unconditioned.direct().underflowed))
                .and_then(|v| v.checked_add(u64::from(pair.conditioned.direct().underflowed)))
                .ok_or(ObservationError::CountOverflow)?;
        }
        let c1 = finite(log_population_count + self.log_mean_frequency)?;
        let c2 = finite(log_add(
            log_population_count + self.log_i00,
            lse(observed_moments.iter().copied())?,
        )?)?;
        Ok(ObjectiveSummary {
            negative_log_likelihood: nll,
            mean_negative_log_likelihood: finite(nll / table.observations as f64)?,
            c1,
            c2,
            log_population_count,
            log_p00: self.log_p00,
            p00: self.p00,
            log_observation_probability: self.log_s,
            raw_prior_log_mass: self.log_mass,
            observations: table.observations,
            distinct_pairs: table.records.len(),
            grid_points: self.grid.frequencies.len(),
            roundoff_corrections: corrections,
            direct_probability_underflows: underflows,
        })
    }
    pub fn evaluate_batch(&self, tables: &[PairTable]) -> Result<Vec<ObjectiveSummary>> {
        let mut output = storage(tables.len())?;
        for table in tables {
            output.push(self.evaluate(table)?);
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn classified_probability_and_allocation_paths() {
        assert_eq!(
            probability(ROUND / 2.0, false).unwrap(),
            LogProbability {
                log_value: 0.0,
                roundoff_corrections: 1
            }
        );
        assert_eq!(
            probability(ROUND * 2.0, false).unwrap_err(),
            ObservationError::InvalidProbability
        );
        for value in [f64::NAN, f64::INFINITY] {
            assert_eq!(
                probability(value, true).unwrap_err(),
                ObservationError::InvalidProbability
            );
        }
        assert_eq!(
            probability(f64::NEG_INFINITY, false).unwrap_err(),
            ObservationError::UnexpectedZeroProbability
        );
        assert_eq!(
            probability(f64::NEG_INFINITY, true).unwrap().direct(),
            DirectProbability {
                value: 0.0,
                underflowed: false
            }
        );
        assert!(probability(-1000.0, false).unwrap().direct().underflowed);
        for value in [0.0, 0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e-20] {
            assert_eq!(
                conditioning(value).unwrap_err(),
                ObservationError::InvalidP00
            );
        }
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.1] {
            assert_eq!(
                valid_denominator(value).unwrap_err(),
                ObservationError::InvalidDenominator
            );
            assert_eq!(
                valid_mass(value).unwrap_err(),
                ObservationError::InconsistentNormalization
            );
        }
        assert_eq!(
            storage::<f64>(usize::MAX).unwrap_err(),
            ObservationError::AllocationFailure
        );
        assert_eq!(
            Grid::log_uniform(0.25, usize::MAX).unwrap_err(),
            ObservationError::AllocationFailure
        );
    }
}
