// SPDX-License-Identifier: GPL-3.0-or-later
//! Conditional posterior kernel for explicitly supplied finite joint priors.
//! Automatic response-prior selection and fitted inference remain scientifically blocked.
use crate::acquisition::AcquisitionError;
use crate::observation::SampleConfiguration;
use std::fmt;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResponseError {
    InvalidPrior,
    DuplicateState,
    NonExchangeable,
    InvalidThreshold,
    InvalidObservation,
    AllImpossible,
    PrecisionLoss,
    InvalidPosterior,
    AllocationFailure,
    NumericalRange,
    Acquisition(AcquisitionError),
}
impl fmt::Display for ResponseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "response error: {self:?}")
    }
}
impl std::error::Error for ResponseError {}
impl From<AcquisitionError> for ResponseError {
    fn from(e: AcquisitionError) -> Self {
        Self::Acquisition(e)
    }
}
type Result<T> = std::result::Result<T, ResponseError>;
fn storage<T>(n: usize) -> Result<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| ResponseError::AllocationFailure)?;
    Ok(v)
}
fn canonical_sum(mut values: Vec<f64>) -> f64 {
    values.sort_unstable_by(f64::total_cmp);
    values.into_iter().sum()
}

#[derive(Debug, Clone, Copy)]
pub struct JointState {
    pub frequency1: f64,
    pub frequency2: f64,
    pub mass: f64,
}
#[derive(Debug, Clone)]
pub struct JointPrior {
    states: Vec<JointState>,
    log_masses: Vec<f64>,
    raw_mass: f64,
}
impl JointPrior {
    pub fn new(states: &[JointState], require_exchangeable: bool) -> Result<Self> {
        if states.is_empty() {
            return Err(ResponseError::InvalidPrior);
        }
        let mut copy = storage(states.len())?;
        let mut masses = storage(states.len())?;
        for (i, &state) in states.iter().enumerate() {
            if !state.frequency1.is_finite()
                || !state.frequency2.is_finite()
                || !(0.0 < state.frequency1
                    && state.frequency1 <= 1.0
                    && 0.0 < state.frequency2
                    && state.frequency2 <= 1.0)
                || !state.mass.is_finite()
                || state.mass <= 0.0
            {
                return Err(ResponseError::InvalidPrior);
            }
            if states[..i]
                .iter()
                .any(|s| s.frequency1 == state.frequency1 && s.frequency2 == state.frequency2)
            {
                return Err(ResponseError::DuplicateState);
            }
            masses.push(state.mass);
            copy.push(state);
        }
        let total = canonical_sum(masses);
        if !total.is_finite() || (total - 1.0).abs() > 1e-12 {
            return Err(ResponseError::InvalidPrior);
        }
        if require_exchangeable {
            for s in states {
                if !states.iter().any(|r| {
                    r.frequency1 == s.frequency2 && r.frequency2 == s.frequency1 && r.mass == s.mass
                }) {
                    return Err(ResponseError::NonExchangeable);
                }
            }
        }
        let mut log_masses = storage(states.len())?;
        for s in states {
            log_masses.push(s.mass.ln() - total.ln());
        }
        Ok(Self {
            states: copy,
            log_masses,
            raw_mass: total,
        })
    }
    pub fn states(&self) -> &[JointState] {
        &self.states
    }
    pub fn raw_mass(&self) -> f64 {
        self.raw_mass
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResponseCall {
    Expanded,
    Contracted,
    NoQualifiedChange,
    Unresolved,
    InvalidModelState,
}
fn qualified_call(regions: [f64; 3], threshold: f64) -> ResponseCall {
    let mut call = ResponseCall::Unresolved;
    let mut qualified = 0;
    for (p, c) in regions.into_iter().zip([
        ResponseCall::Expanded,
        ResponseCall::Contracted,
        ResponseCall::NoQualifiedChange,
    ]) {
        if p >= threshold {
            qualified += 1;
            call = c;
        }
    }
    if qualified == 1 {
        call
    } else {
        ResponseCall::Unresolved
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Posterior {
    pub expanded: f64,
    pub contracted: f64,
    pub no_change: f64,
    pub log_evidence: f64,
    pub mean_log_response: f64,
    pub credible_low: f64,
    pub credible_high: f64,
    pub state_probabilities: Vec<f64>,
    pub state_log_probabilities: Vec<f64>,
    pub posterior_underflows: usize,
    pub call: ResponseCall,
    pub raw_prior_mass: f64,
    pub posterior_mass: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseRecord {
    pub count1: u64,
    pub count2: u64,
    pub result: Result<Posterior>,
}
impl ResponseRecord {
    pub fn call(&self) -> ResponseCall {
        match &self.result {
            Ok(p) => p.call,
            Err(_) => ResponseCall::InvalidModelState,
        }
    }
}
#[derive(Debug, Clone)]
pub struct ResponseModel {
    prior: JointPrior,
    sample1: SampleConfiguration,
    sample2: SampleConfiguration,
    delta: f64,
    threshold: f64,
    responses: Vec<f64>,
    response_order: Vec<usize>,
    frequency_indices: Vec<[usize; 2]>,
    frequency_counts: [usize; 2],
}
impl ResponseModel {
    pub fn new(
        prior: JointPrior,
        sample1: SampleConfiguration,
        sample2: SampleConfiguration,
        delta: f64,
        threshold: f64,
    ) -> Result<Self> {
        if !delta.is_finite()
            || delta < 0.0
            || !threshold.is_finite()
            || threshold <= 0.5
            || threshold > 1.0
        {
            return Err(ResponseError::InvalidThreshold);
        }
        for s in &prior.states {
            sample1
                .acquisition
                .log_pmf(0, s.frequency1, sample1.depth)?;
            sample2
                .acquisition
                .log_pmf(0, s.frequency2, sample2.depth)?;
        }
        // Exact binary64 frequency identity. Build after original validation
        // so allocation cannot mask any pre-existing acquisition error.
        let mut unique: [Vec<u64>; 2] =
            [storage(prior.states.len())?, storage(prior.states.len())?];
        for state in &prior.states {
            unique[0].push(state.frequency1.to_bits());
            unique[1].push(state.frequency2.to_bits());
        }
        for values in &mut unique {
            values.sort_unstable();
            values.dedup();
        }
        let frequency_counts = [unique[0].len(), unique[1].len()];
        let mut frequency_indices = storage(prior.states.len())?;
        for state in &prior.states {
            frequency_indices.push([
                unique[0]
                    .binary_search(&state.frequency1.to_bits())
                    .expect("inserted frequency"),
                unique[1]
                    .binary_search(&state.frequency2.to_bits())
                    .expect("inserted frequency"),
            ]);
        }
        let mut responses = storage(prior.states.len())?;
        let mut response_order = storage(prior.states.len())?;
        for (index, state) in prior.states.iter().enumerate() {
            responses.push(state.frequency2.ln() - state.frequency1.ln());
            response_order.push(index);
        }
        response_order
            .sort_unstable_by(|&a, &b| responses[a].total_cmp(&responses[b]).then(a.cmp(&b)));
        Ok(Self {
            prior,
            sample1,
            sample2,
            delta,
            threshold,
            responses,
            response_order,
            frequency_indices,
            frequency_counts,
        })
    }
    pub fn infer(&self, count1: u64, count2: u64) -> Result<Posterior> {
        self.infer_masked(count1, count2, None)
    }
    /// Opt-in omission certification. The accepted zero-omission calculation
    /// is returned without changing any arithmetic or posterior field.
    pub fn infer_certified(
        &self,
        count1: u64,
        count2: u64,
    ) -> Result<(Posterior, crate::certification::Certificate)> {
        self.infer_certified_detailed(count1, count2)
            .map_err(|e| e.error)
    }
    /// Detailed failure classification; successful calculations are unchanged.
    pub fn infer_certified_detailed(
        &self,
        count1: u64,
        count2: u64,
    ) -> std::result::Result<
        (Posterior, crate::certification::Certificate),
        crate::certification::CertificationFailure,
    > {
        match self.infer(count1, count2) {
            Ok(p) => {
                return Ok((
                    p,
                    crate::certification::Certificate::exact(self.prior.states.len()),
                ))
            }
            Err(ResponseError::PrecisionLoss) => {}
            Err(e) => return Err(e.into()),
        }
        let samples = [self.sample1, self.sample2];
        let counts = [count1, count2];
        let (mask, mut certificate) =
            crate::certification::certify(&self.prior.states, samples, counts)?;
        let posterior = self.infer_masked(count1, count2, Some(&mask))?;
        crate::certification::certify_intervals(
            &self.prior.states,
            samples,
            counts,
            &mask,
            &posterior,
            &mut certificate,
        )?;
        Ok((posterior, certificate))
    }
    fn infer_masked(
        &self,
        count1: u64,
        count2: u64,
        omissions: Option<&[bool]>,
    ) -> Result<Posterior> {
        if count1 == 0 && count2 == 0 {
            return Err(ResponseError::InvalidObservation);
        }
        let mut likelihoods = storage(self.prior.states.len())?;
        let mut cache1 = storage(self.frequency_counts[0])?;
        let mut cache2 = storage(self.frequency_counts[1])?;
        cache1.resize(self.frequency_counts[0], None);
        cache2.resize(self.frequency_counts[1], None);
        let mut max = f64::NEG_INFINITY;
        for (i, s) in self.prior.states.iter().enumerate() {
            if omissions.is_some_and(|mask| mask[i]) {
                likelihoods.push(f64::NEG_INFINITY);
                continue;
            }
            // Lazy first retained occurrence preserves original state/axis
            // acquisition error order. Keep BOTH calls before the precision
            // guard, exactly as in the original kernel.
            let [i1, i2] = self.frequency_indices[i];
            let a = match cache1[i1] {
                Some(value) => value,
                None => {
                    let value = self.sample1.acquisition.log_pmf(
                        count1,
                        s.frequency1,
                        self.sample1.depth,
                    )?;
                    cache1[i1] = Some(value);
                    value
                }
            };
            let b = match cache2[i2] {
                Some(value) => value,
                None => {
                    let value = self.sample2.acquisition.log_pmf(
                        count2,
                        s.frequency2,
                        self.sample2.depth,
                    )?;
                    cache2[i2] = Some(value);
                    value
                }
            };
            for log in [a, b] {
                if log.is_finite() && 16.0 * f64::EPSILON * log.abs() > 5e-10 {
                    return Err(ResponseError::PrecisionLoss);
                }
            }
            let log = a + b;
            if log.is_nan() || log == f64::INFINITY {
                return Err(ResponseError::NumericalRange);
            }
            likelihoods.push(log);
            max = max.max(log);
        }
        if max == f64::NEG_INFINITY {
            return Err(ResponseError::AllImpossible);
        }
        // Preserve prior weights before reconstructing the potentially large log evidence.
        let mut weights = storage(likelihoods.len())?;
        let mut offset = f64::NEG_INFINITY;
        for (l, p) in likelihoods.iter().zip(&self.prior.log_masses) {
            let v = (l - max) + p;
            offset = offset.max(v);
            weights.push(v);
        }
        let mut scaled = storage(weights.len())?;
        scaled.extend(weights.iter().map(|v| (v - offset).exp()));
        let norm = canonical_sum(scaled);
        if !norm.is_finite() || norm <= 0.0 {
            return Err(ResponseError::InvalidPosterior);
        }
        let mut probabilities = storage(weights.len())?;
        let mut log_probabilities = storage(weights.len())?;
        let mut underflows = 0;
        let mut regions = [
            storage(weights.len())?,
            storage(weights.len())?,
            storage(weights.len())?,
        ];
        let mut mean = 0.0;
        for (i, &log) in weights.iter().enumerate() {
            let log_p = (log - offset) - norm.ln();
            let p = log_p.exp();
            underflows += usize::from(p == 0.0 && log_p.is_finite());
            log_probabilities.push(log_p);
            let response = self.responses[i];
            if !p.is_finite() || !(0.0..=1.0).contains(&p) {
                return Err(ResponseError::InvalidPosterior);
            }
            if response > self.delta {
                regions[0].push(p);
            } else if response < -self.delta {
                regions[1].push(p);
            } else {
                regions[2].push(p);
            }
            mean += p * response;
            probabilities.push(p);
        }
        let [exp_terms, con_terms, no_terms] = regions;
        let expanded = canonical_sum(exp_terms);
        let contracted = canonical_sum(con_terms);
        let no_change = canonical_sum(no_terms);
        let total = expanded + contracted + no_change;
        if !total.is_finite()
            || (total - 1.0).abs() > 3e-12
            || [expanded, contracted, no_change]
                .iter()
                .any(|p| *p < 0.0 || *p > 1.0 + 3e-12)
        {
            return Err(ResponseError::InvalidPosterior);
        }
        let quantile = |q: f64| {
            let mut cumulative = 0.0;
            for &index in &self.response_order {
                cumulative += probabilities[index];
                if cumulative >= q {
                    return self.responses[index];
                }
            }
            self.responses[*self
                .response_order
                .last()
                .expect("validated nonempty prior")]
        };
        let call = qualified_call([expanded, contracted, no_change], self.threshold);
        let log_evidence = max + offset + norm.ln();
        if !log_evidence.is_finite() || log_evidence > 3e-12 || !mean.is_finite() {
            return Err(ResponseError::NumericalRange);
        }
        Ok(Posterior {
            expanded,
            contracted,
            no_change,
            log_evidence,
            mean_log_response: mean,
            credible_low: quantile(0.025),
            credible_high: quantile(0.975),
            state_probabilities: probabilities,
            state_log_probabilities: log_probabilities,
            posterior_underflows: underflows,
            call,
            raw_prior_mass: self.prior.raw_mass,
            posterior_mass: total,
        })
    }
    pub fn assess(&self, count1: u64, count2: u64) -> ResponseRecord {
        ResponseRecord {
            count1,
            count2,
            result: self.infer(count1, count2),
        }
    }
    pub fn assess_batch(&self, counts: &[(u64, u64)]) -> Result<Vec<ResponseRecord>> {
        let mut results = storage(counts.len())?;
        for &(a, b) in counts {
            results.push(self.assess(a, b));
        }
        Ok(results)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn allocation_failure_is_classified() {
        let threshold = f64::from_bits(0.5_f64.to_bits() + 1);
        assert_eq!(
            qualified_call([threshold, threshold, 0.0], threshold),
            ResponseCall::Unresolved
        );
        assert_eq!(
            storage::<f64>(usize::MAX).unwrap_err(),
            ResponseError::AllocationFailure
        );
    }
}
