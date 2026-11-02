// SPDX-License-Identifier: GPL-3.0-or-later
//! Rigorous cumulative Poisson omission bounds and outward CDF enclosures.
//! Raw prior normalization cancels from the evidence ratio. No guarded log-PMF
//! is used to prove negligibility or the retained evidence lower bound.
use crate::acquisition::Acquisition;
use crate::observation::SampleConfiguration;
use crate::response::{JointState, Posterior, ResponseError};
use serde::Serialize;

pub const MASS_EXPONENT: i64 = -80;
pub const MEAN_EXPONENT: i64 = -64;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum CertificationReason {
    RetainedEvidenceAnchorUnavailable,
    CumulativeOmissionBoundTooLoose,
    UnsupportedArithmeticDomain,
    AmbiguousIntervalCertification,
    InvalidCertificationConfiguration,
    InternalCertificationFailure,
}
impl CertificationReason {
    pub fn message(self) -> &'static str {
        match self {
            Self::RetainedEvidenceAnchorUnavailable => {
                "No retained state supplies a certified evidence lower bound; treat this result as unresolved."
            }
            Self::CumulativeOmissionBoundTooLoose => {
                "The cumulative omission bound exceeds the fixed budget; treat this result as unresolved."
            }
            Self::UnsupportedArithmeticDomain => {
                "Certificate arithmetic is outside its supported domain; inspect inputs and supported numeric ranges."
            }
            Self::AmbiguousIntervalCertification => {
                "A credible interval endpoint could not be certified; treat this result as unresolved."
            }
            Self::InvalidCertificationConfiguration => {
                "The certification configuration is unsupported; check the model and count limits."
            }
            Self::InternalCertificationFailure => "An internal certification invariant failed; retain inputs and report the failure.",
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CertificationFailure {
    pub error: ResponseError,
    pub reason: Option<CertificationReason>,
}
impl CertificationFailure {
    fn precision(reason: CertificationReason) -> Self {
        Self {
            error: ResponseError::PrecisionLoss,
            reason: Some(reason),
        }
    }
}
impl From<ResponseError> for CertificationFailure {
    fn from(error: ResponseError) -> Self {
        Self {
            error,
            reason: (error == ResponseError::PrecisionLoss)
                .then_some(CertificationReason::InternalCertificationFailure),
        }
    }
}
impl From<crate::acquisition::AcquisitionError> for CertificationFailure {
    fn from(error: crate::acquisition::AcquisitionError) -> Self {
        ResponseError::from(error).into()
    }
}
type Result<T> = std::result::Result<T, CertificationFailure>;
fn fail<T>() -> Result<T> {
    Err(CertificationFailure::precision(
        CertificationReason::UnsupportedArithmeticDomain,
    ))
}
fn reason<T>(reason: CertificationReason) -> Result<T> {
    Err(CertificationFailure::precision(reason))
}
fn buffer<T>(capacity: usize) -> Result<Vec<T>> {
    let mut value = Vec::new();
    value
        .try_reserve_exact(capacity)
        .map_err(|_| CertificationFailure::from(ResponseError::AllocationFailure))?;
    Ok(value)
}
fn collect<T>(values: impl Iterator<Item = T>, capacity: usize) -> Result<Vec<T>> {
    let mut result = buffer(capacity)?;
    result.extend(values);
    Ok(result)
}
fn down(x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else {
        f64::from_bits(if x > 0.0 {
            x.to_bits() - 1
        } else {
            x.to_bits() + 1
        })
    }
}
fn up(x: f64) -> f64 {
    if x == 0.0 {
        f64::from_bits(1)
    } else if x == f64::INFINITY {
        x
    } else {
        f64::from_bits(if x > 0.0 {
            x.to_bits() + 1
        } else {
            x.to_bits() - 1
        })
    }
}
fn floor2(x: f64) -> i64 {
    let bits = x.to_bits();
    let e = ((bits >> 52) & 2047) as i64;
    if e == 0 {
        63 - bits.leading_zeros() as i64 - 1074
    } else {
        e - 1023
    }
}
fn ceil2_integer(n: u64) -> i64 {
    if n <= 1 {
        0
    } else {
        64 - (n - 1).leading_zeros() as i64
    }
}
fn mass_exponents(x: f64) -> (i64, i64) {
    let lo = floor2(x);
    let power = if x.is_normal() {
        x.to_bits() & ((1_u64 << 52) - 1) == 0
    } else {
        x.to_bits().count_ones() == 1
    };
    (lo, lo + i64::from(!power))
}
fn mean_bounds(depth: f64, frequency: f64) -> Result<(f64, f64)> {
    let mu = depth * frequency;
    if !mu.is_finite() || mu < 0.0 {
        return fail();
    }
    let lo = down(mu);
    let hi = up(mu);
    if !hi.is_finite() {
        return fail();
    }
    Ok((lo, hi))
}
fn pmf_upper(n: u64, mu_lo: f64) -> Result<i64> {
    if mu_lo > (i64::MAX / 4) as f64 {
        return fail();
    }
    Ok((n as i64 - (mu_lo * 0.5).floor() as i64).min(0))
}
fn pmf_lower(n: u64, mu_lo: f64, mu_hi: f64) -> Option<i64> {
    if mu_hi > (i64::MAX / 4) as f64 {
        return None;
    }
    let base = -(2.0 * mu_hi).ceil() as i64;
    if n == 0 {
        return Some(base);
    }
    if mu_lo <= 0.0 {
        return None;
    }
    base.checked_add((n as i64).checked_mul(floor2(mu_lo) - ceil2_integer(n))?)
}
fn cumulative_bound(maximum: i64, count: usize, anchor: i64) -> Result<(i64, i64)> {
    let cumulative = maximum
        .checked_add(ceil2_integer(count as u64))
        .ok_or_else(|| {
            CertificationFailure::precision(CertificationReason::UnsupportedArithmeticDomain)
        })?;
    let ratio = cumulative.checked_sub(anchor).ok_or_else(|| {
        CertificationFailure::precision(CertificationReason::UnsupportedArithmeticDomain)
    })?;
    if ratio > MASS_EXPONENT || ratio + 11 > MEAN_EXPONENT {
        return reason(CertificationReason::CumulativeOmissionBoundTooLoose);
    }
    Ok((cumulative, ratio))
}

#[derive(Debug, Clone, Serialize)]
pub struct Certificate {
    pub method: &'static str,
    pub total_states: usize,
    pub retained_states: usize,
    pub omitted_states: usize,
    pub omitted_indices: Vec<usize>,
    pub maximum_omitted_raw_weight_log2_upper: Option<i64>,
    pub maximum_omitted_log_contribution_upper: Option<f64>,
    pub cumulative_omitted_raw_weight_log2_upper: Option<i64>,
    pub retained_evidence_raw_log2_lower: Option<i64>,
    pub posterior_omitted_mass_log2_upper: Option<i64>,
    pub certified_posterior_mass_upper: f64,
    pub certified_region_error_upper: f64,
    pub certified_relative_and_log_evidence_error_upper: f64,
    pub certified_mean_error_upper: f64,
    pub requested_mass_log2_budget: i64,
    pub requested_mean_log2_budget: i64,
    pub quantile_cdf_margins: Option<[[f64; 2]; 2]>,
    pub retained_roundoff: &'static str,
}
impl Certificate {
    pub fn exact(n: usize) -> Self {
        Self {
            method: "dyadic Poisson pgf/anchor; outward scaled CDF",
            total_states: n,
            retained_states: n,
            omitted_states: 0,
            omitted_indices: vec![],
            maximum_omitted_raw_weight_log2_upper: None,
            maximum_omitted_log_contribution_upper: None,
            cumulative_omitted_raw_weight_log2_upper: None,
            retained_evidence_raw_log2_lower: None,
            posterior_omitted_mass_log2_upper: None,
            certified_posterior_mass_upper: 0.0,
            certified_region_error_upper: 0.0,
            certified_relative_and_log_evidence_error_upper: 0.0,
            certified_mean_error_upper: 0.0,
            requested_mass_log2_budget: MASS_EXPONENT,
            requested_mean_log2_budget: MEAN_EXPONENT,
            quantile_cdf_margins: None,
            retained_roundoff:
                "unchanged accepted calculation; guard is not a universal accuracy theorem",
        }
    }
}

pub(crate) fn certify(
    states: &[JointState],
    samples: [SampleConfiguration; 2],
    counts: [u64; 2],
) -> Result<(Vec<bool>, Certificate)> {
    if samples
        .iter()
        .any(|s| !matches!(s.acquisition, Acquisition::Poisson))
        || counts.iter().any(|n| *n > 1_000_000)
    {
        return reason(CertificationReason::InvalidCertificationConfiguration);
    }
    let mut mask = buffer(states.len())?;
    let mut omitted = buffer(states.len())?;
    let mut maximum = i64::MIN;
    let mut anchor = i64::MIN;
    for (i, s) in states.iter().enumerate() {
        let response = s.frequency2.ln() - s.frequency1.ln();
        if !response.is_finite() || response.abs() > 1024.0 {
            return reason(CertificationReason::InvalidCertificationConfiguration);
        }
        let f = [s.frequency1, s.frequency2];
        let mut unsafe_state = false;
        let (mut lower, mut upper) = mass_exponents(s.mass);
        let mut lower_possible = true;
        for j in 0..2 {
            let log = samples[j]
                .acquisition
                .log_pmf(counts[j], f[j], samples[j].depth)?;
            if log.is_nan() || log == f64::INFINITY {
                return fail();
            }
            unsafe_state |= log.is_finite() && 16.0 * f64::EPSILON * log.abs() > 5e-10;
            let (lo, hi) = mean_bounds(samples[j].depth, f[j])?;
            upper = upper
                .checked_add(pmf_upper(counts[j], lo)?)
                .ok_or_else(|| {
                    CertificationFailure::precision(
                        CertificationReason::UnsupportedArithmeticDomain,
                    )
                })?;
            if let Some(l) = pmf_lower(counts[j], lo, hi) {
                lower = lower.checked_add(l).ok_or_else(|| {
                    CertificationFailure::precision(
                        CertificationReason::UnsupportedArithmeticDomain,
                    )
                })?;
            } else {
                lower_possible = false;
            }
        }
        mask.push(unsafe_state);
        if unsafe_state {
            omitted.push(i);
            maximum = maximum.max(upper);
        } else if lower_possible {
            anchor = anchor.max(lower);
        }
    }
    if omitted.is_empty() {
        return Ok((mask, Certificate::exact(states.len())));
    }
    if anchor == i64::MIN {
        return reason(CertificationReason::RetainedEvidenceAnchorUnavailable);
    }
    // Keep every previously accepted coarse certificate byte-for-byte. A
    // failed coarse anchor can be sharpened by the existing outward scaled
    // arithmetic; no guarded PMF logarithm proves this lower bound.
    let mut refined_used = false;
    if matches!(
        cumulative_bound(maximum, omitted.len(), anchor),
        Err(CertificationFailure {
            reason: Some(CertificationReason::CumulativeOmissionBoundTooLoose),
            ..
        })
    ) {
        if let Some((refined_maximum, refined_anchor)) =
            refined_bounds(states, samples, counts, &mask)
        {
            maximum = maximum.min(refined_maximum);
            anchor = anchor.max(refined_anchor);
            refined_used = true;
        }
    }
    let (cumulative, ratio) = cumulative_bound(maximum, omitted.len(), anchor)?;
    // A positive, conservative representable upper bound even for exponents
    // far below the floating range. Never use exp underflow as a certificate.
    let display = 2_f64.powi(ratio.max(-1022) as i32);
    let mut c = Certificate::exact(states.len());
    if refined_used {
        c.method = "dyadic Poisson scaled PMF refinement; outward scaled CDF";
    }
    c.retained_states -= omitted.len();
    c.omitted_states = omitted.len();
    c.omitted_indices = omitted;
    c.maximum_omitted_raw_weight_log2_upper = Some(maximum);
    // ln(2)<3/4; for negative exponent multiplying by 1/2 is an upper
    // natural-log envelope (ln(2)>1/2). Integer dyadic bound is authoritative.
    c.maximum_omitted_log_contribution_upper = Some(if maximum < 0 {
        maximum as f64 * 0.5
    } else {
        maximum as f64
    });
    c.cumulative_omitted_raw_weight_log2_upper = Some(cumulative);
    c.retained_evidence_raw_log2_lower = Some(anchor);
    c.posterior_omitted_mass_log2_upper = Some(ratio);
    c.certified_posterior_mass_upper = display;
    c.certified_region_error_upper = display;
    c.certified_relative_and_log_evidence_error_upper = display;
    c.certified_mean_error_upper = display * 2048.0;
    Ok((mask, c))
}

/// Positive outward interval with an unbounded (i64) binary scale. Every
/// arithmetic operation is enclosed by adjacent IEEE754 binary64 numbers.
#[derive(Debug, Clone, Copy)]
struct Scaled {
    lo: f64,
    hi: f64,
    exponent: i64,
}
impl Scaled {
    fn new(lo: f64, hi: f64, exponent: i64) -> Result<Self> {
        if !(lo >= 0.0 && hi > 0.0 && lo <= hi && hi.is_finite()) {
            return fail();
        }
        let shift = floor2(hi);
        let factor = 2_f64.powi(-shift as i32);
        if !factor.is_finite() || factor == 0.0 {
            return fail();
        }
        Ok(Self {
            lo: down(lo * factor),
            hi: up(hi * factor),
            exponent: exponent.checked_add(shift).ok_or_else(|| {
                CertificationFailure::precision(CertificationReason::UnsupportedArithmeticDomain)
            })?,
        })
    }
    fn from(x: f64) -> Result<Self> {
        if !(x > 0.0 && x.is_finite()) {
            return fail();
        }
        let bits = x.to_bits();
        let e = ((bits >> 52) & 2047) as i64;
        if e != 0 {
            let m = f64::from_bits((1023_u64 << 52) | (bits & ((1_u64 << 52) - 1)));
            Ok(Self {
                lo: m,
                hi: m,
                exponent: e - 1023,
            })
        } else {
            let shift = 63 - bits.leading_zeros() as i64;
            let m = bits as f64 / 2_f64.powi(shift as i32);
            Ok(Self {
                lo: m,
                hi: m,
                exponent: shift - 1074,
            })
        }
    }
    fn mul(self, b: Self) -> Result<Self> {
        Self::new(
            down(self.lo * b.lo),
            up(self.hi * b.hi),
            self.exponent.checked_add(b.exponent).ok_or_else(|| {
                CertificationFailure::precision(CertificationReason::UnsupportedArithmeticDomain)
            })?,
        )
    }
    fn power(self, mut n: u64) -> Result<Self> {
        let mut base = self;
        let mut result = Self::from(1.0)?;
        while n > 0 {
            if n & 1 != 0 {
                result = result.mul(base)?;
            }
            n >>= 1;
            if n > 0 {
                base = base.mul(base)?;
            }
        }
        Ok(result)
    }
    fn aligned(self, e: i64) -> (f64, f64) {
        let gap = e - self.exponent;
        if gap >= 1023 {
            (0.0, f64::MIN_POSITIVE)
        } else {
            let factor = 2_f64.powi(-(gap as i32));
            (down(self.lo * factor), up(self.hi * factor))
        }
    }
    fn add(self, b: Self) -> Result<Self> {
        let e = self.exponent.max(b.exponent);
        let a = self.aligned(e);
        let b = b.aligned(e);
        Self::new(down(a.0 + b.0), up(a.1 + b.1), e)
    }
    fn ratio(self, b: Self) -> Result<(f64, f64)> {
        if b.lo <= 0.0 {
            return fail();
        }
        let gap = self.exponent - b.exponent;
        if gap < -1022 {
            let bound = up(up(self.hi / b.lo) * f64::MIN_POSITIVE);
            if !bound.is_finite() {
                return fail();
            }
            return Ok((0.0, bound));
        }
        if gap > 1022 {
            return fail();
        }
        let factor = 2_f64.powi(gap as i32);
        let lo = down(down(self.lo / b.hi) * factor);
        let hi = up(up(self.hi / b.lo) * factor);
        if !lo.is_finite() || !hi.is_finite() {
            return fail();
        }
        Ok((lo, hi))
    }
}
/// Outward division retains the unbounded binary scale. All denominators here
/// are strictly positive factorial intervals, never observation normalizers.
fn scaled_divide(a: Scaled, b: Scaled) -> Result<Scaled> {
    if b.lo <= 0.0 {
        return fail();
    }
    Scaled::new(
        down(a.lo / b.hi),
        up(a.hi / b.lo),
        a.exponent.checked_sub(b.exponent).ok_or_else(|| {
            CertificationFailure::precision(CertificationReason::UnsupportedArithmeticDomain)
        })?,
    )
}
fn factorial_interval(n: u64) -> Result<Scaled> {
    let mut result = Scaled::from(1.0)?;
    for k in 2..=n {
        // certify() restricts n to <=1,000,000, hence this conversion is exact.
        result = result.mul(Scaled::from(k as f64)?)?;
    }
    Ok(result)
}
fn refined_pmf_interval(n: u64, lo: f64, hi: f64, factorial: Scaled) -> Option<Scaled> {
    let exponential = exp_negative(lo, hi).ok()?;
    let power = if n == 0 {
        Scaled::from(1.0).ok()?
    } else {
        if lo <= 0.0 {
            return None;
        }
        Scaled::new(lo, hi, 0).ok()?.power(n).ok()?
    };
    scaled_divide(exponential.mul(power).ok()?, factorial).ok()
}
fn refined_bounds(
    states: &[JointState],
    samples: [SampleConfiguration; 2],
    counts: [u64; 2],
    mask: &[bool],
) -> Option<(i64, i64)> {
    let factorials = [
        factorial_interval(counts[0]).ok()?,
        factorial_interval(counts[1]).ok()?,
    ];
    let (mut maximum, mut anchor) = (i64::MIN, i64::MIN);
    for (state, omitted) in states.iter().zip(mask) {
        let frequencies = [state.frequency1, state.frequency2];
        let (mut lower, mut upper) = mass_exponents(state.mass);
        let mut lower_possible = true;
        for j in 0..2 {
            let (lo, hi) = mean_bounds(samples[j].depth, frequencies[j]).ok()?;
            if let Some(p) = refined_pmf_interval(counts[j], lo, hi, factorials[j]) {
                upper = upper.checked_add(p.exponent.checked_add(mass_exponents(p.hi).1)?)?;
                if p.lo > 0.0 {
                    lower = lower.checked_add(p.exponent.checked_add(floor2(p.lo))?)?;
                } else {
                    lower_possible = false;
                }
            } else {
                // An unavailable refinement retains the original rigorous bound.
                upper = upper.checked_add(pmf_upper(counts[j], lo).ok()?)?;
                if let Some(l) = pmf_lower(counts[j], lo, hi) {
                    lower = lower.checked_add(l)?;
                } else {
                    lower_possible = false;
                }
            }
        }
        if *omitted {
            maximum = maximum.max(upper);
        } else if lower_possible {
            anchor = anchor.max(lower);
        }
    }
    if maximum == i64::MIN || anchor == i64::MIN {
        None
    } else {
        Some((maximum, anchor))
    }
}

fn sum(mut x: Vec<Scaled>) -> Result<Option<Scaled>> {
    x.sort_unstable_by(|a, b| {
        a.exponent
            .cmp(&b.exponent)
            .then(a.lo.total_cmp(&b.lo))
            .then(a.hi.total_cmp(&b.hi))
    });
    let mut result = None;
    for s in x {
        result = Some(match result {
            None => s,
            Some(a) => Scaled::add(a, s)?,
        });
    }
    Ok(result)
}
fn exp_negative(mut lo: f64, mut hi: f64) -> Result<Scaled> {
    let mut squares = 0;
    while hi > 0.5 {
        lo *= 0.5;
        hi *= 0.5;
        squares += 1;
        if squares > 1024 {
            return fail();
        }
    }
    let (mut tl, mut th, mut sl, mut sh) = (1.0, 1.0, 1.0, 1.0);
    for k in 1..=24 {
        tl = down(down(tl * lo) / k as f64);
        th = up(up(th * hi) / k as f64);
        sl = down(sl + tl);
        sh = up(sh + th);
    }
    let next = up(up(th * hi) / 25.0);
    let remainder = up(next / down(1.0 - up(hi / 26.0)));
    let mut value = Scaled::new(down(1.0 / up(sh + remainder)), up(1.0 / sl), 0)?;
    for _ in 0..squares {
        value = value.mul(value)?;
    }
    Ok(value)
}

pub(crate) fn certify_intervals(
    states: &[JointState],
    samples: [SampleConfiguration; 2],
    counts: [u64; 2],
    mask: &[bool],
    posterior: &Posterior,
    certificate: &mut Certificate,
) -> Result<()> {
    if certificate.omitted_states == 0 {
        return Ok(());
    }
    let mut weights = buffer(states.len())?;
    // Identical per-axis frequencies have identical interval arithmetic. Keep
    // first-encounter evaluation and all subsequent joint sums unchanged.
    // Capacity bounds keep every new allocation fallible and classified.
    let mut factors_cache: [Vec<(u64, Scaled)>; 2] = [buffer(states.len())?, buffer(states.len())?];
    for (i, s) in states.iter().enumerate() {
        if mask[i] {
            continue;
        }
        let f = [s.frequency1, s.frequency2];
        let mut factors = [None, None];
        for j in 0..2 {
            let weight = match factors_cache[j].binary_search_by_key(&f[j].to_bits(), |x| x.0) {
                Ok(index) => factors_cache[j][index].1,
                Err(index) => {
                    let (lo, hi) = mean_bounds(samples[j].depth, f[j])?;
                    let weight = Scaled::from(f[j])?
                        .power(counts[j])?
                        .mul(exp_negative(lo, hi)?)?;
                    factors_cache[j].insert(index, (f[j].to_bits(), weight));
                    weight
                }
            };
            factors[j] = Some(weight);
        }
        let weight = factors[0]
            .ok_or_else(|| {
                CertificationFailure::precision(CertificationReason::InternalCertificationFailure)
            })?
            .mul(factors[1].ok_or_else(|| {
                CertificationFailure::precision(CertificationReason::InternalCertificationFailure)
            })?)?
            .mul(Scaled::from(s.mass)?)?;
        weights.push((s.frequency2.ln() - s.frequency1.ln(), weight));
    }
    let total = sum(collect(weights.iter().map(|x| x.1), weights.len())?)?.ok_or_else(|| {
        CertificationFailure::precision(CertificationReason::InternalCertificationFailure)
    })?;
    let rho = certificate.certified_posterior_mass_upper;
    let mut margins = [[0.0; 2]; 2];
    for (i, (q, endpoint)) in [
        (0.025, posterior.credible_low),
        (0.975, posterior.credible_high),
    ]
    .into_iter()
    .enumerate()
    {
        let before = sum(collect(
            weights.iter().filter(|x| x.0 < endpoint).map(|x| x.1),
            weights.len(),
        )?)?;
        let at = sum(collect(
            weights.iter().filter(|x| x.0 <= endpoint).map(|x| x.1),
            weights.len(),
        )?)?
        .ok_or_else(|| {
            CertificationFailure::precision(CertificationReason::InternalCertificationFailure)
        })?;
        let upper_before = match before {
            None => 0.0,
            Some(s) => s.ratio(total)?.1,
        };
        let lower_at = at.ratio(total)?.0;
        if up(upper_before + rho) >= q || down(lower_at - rho) < q {
            return reason(CertificationReason::AmbiguousIntervalCertification);
        }
        margins[i] = [
            down(q - up(upper_before + rho)),
            down(down(lower_at - rho) - q),
        ];
    }
    certificate.quantile_cdf_margins = Some(margins);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn arithmetic_domain_and_internal_invariant_have_distinct_reasons() {
        assert_eq!(
            mean_bounds(f64::MAX, 2.0).unwrap_err().reason,
            Some(CertificationReason::UnsupportedArithmeticDomain)
        );
        let states = [JointState {
            frequency1: 0.1,
            frequency2: 0.1,
            mass: 1.0,
        }];
        let samples = [SampleConfiguration {
            acquisition: Acquisition::Poisson,
            depth: 1.0,
        }; 2];
        let prior = crate::response::JointPrior::new(&states, true).unwrap();
        let model =
            crate::response::ResponseModel::new(prior, samples[0], samples[1], 0.3, 1.0).unwrap();
        let posterior = model.infer(1, 1).unwrap();
        let mut certificate = Certificate::exact(1);
        certificate.omitted_states = 1;
        assert_eq!(
            certify_intervals(
                &states,
                samples,
                [1, 1],
                &[true],
                &posterior,
                &mut certificate
            )
            .unwrap_err()
            .reason,
            Some(CertificationReason::InternalCertificationFailure)
        );
    }
    #[test]
    fn cumulative_many_tiny_states_are_not_individually_accepted() {
        assert!(cumulative_bound(-81, 1, 0).is_ok());
        assert!(cumulative_bound(-81, 4, 0).is_err());
        assert!(cumulative_bound(-83, 4, 0).is_ok());
        assert!(cumulative_bound(-80, 2, 0).is_err());
    }
    #[test]
    fn binary_envelopes_include_subnormal_prior_masses() {
        assert_eq!(mass_exponents(f64::from_bits(1)), (-1074, -1074));
        assert_eq!(mass_exponents(f64::from_bits(3)), (-1073, -1072));
        assert_eq!(mass_exponents(0.5), (-1, -1));
        assert_eq!(mass_exponents(0.75), (-1, 0));
    }
    #[test]
    fn scaled_intervals_do_not_underflow_evidence_or_ignore_small_summands() {
        let a = Scaled::from(1e-6).unwrap().power(3000).unwrap();
        assert!(a.exponent < -50_000 && a.hi > 0.0);
        let b = Scaled::new(1.0, 1.0, a.exponent - 4000).unwrap();
        let total = sum(vec![a, b]).unwrap().unwrap();
        let (lo, hi) = a.ratio(total).unwrap();
        assert!(lo < 1.0 && hi >= 1.0);
        let e = exp_negative(258529.0, 258529.0).unwrap();
        assert!(e.exponent < -300_000 && e.hi > 0.0);
    }
    #[test]
    fn nonfinite_interval_is_refused() {
        assert!(mean_bounds(f64::MAX, 2.0).is_err());
        assert!(Scaled::from(f64::INFINITY).is_err());
        assert!(Scaled::from(0.0).is_err());
    }
    #[test]
    fn subnormal_ratio_encloses_a_denominator_lower_mantissa_below_one() {
        let a = Scaled::new(1.5, 1.5, -1023).unwrap();
        let b = Scaled::new(0.5, 1.0, 0).unwrap();
        let (lo, hi) = a.ratio(b).unwrap();
        let exact_upper = 3.0 * 2_f64.powi(-1023);
        assert_eq!(lo, 0.0);
        assert!(hi >= exact_upper);
    }
    #[test]
    fn wide_interval_rescaling_does_not_round_a_lower_endpoint_up() {
        let s = Scaled::new(f64::from_bits(3), 4.0, 0).unwrap();
        assert_eq!(s.lo, 0.0);
        assert!(s.hi >= 1.0);
    }
}

#[cfg(test)]
mod refined_tests {
    use super::*;
    #[test]
    fn high_count_scaled_pmf_matches_independent_decimal_dyadic_envelopes() {
        // Independent110-digit Decimal factorial intervals and log probabilities.
        for (n, mu, exponent) in [
            (90_000, 90_000.0, -10),
            (90_000, 0.5, -1_441_356),
            (90_000, 180_000.0, -39_853),
            (1_000_000, 1_000_000.0, -12),
        ] {
            let (lo, hi) = mean_bounds(mu, 1.0).unwrap();
            let p = refined_pmf_interval(n, lo, hi, factorial_interval(n).unwrap()).unwrap();
            assert!(p.lo > 0.0 && p.lo <= p.hi);
            assert_eq!(p.exponent, exponent);
            assert!(p.hi / p.lo < 1.000_001);
        }
    }
    #[test]
    fn high_count_low_and_high_mean_omissions_get_rigorous_joint_certificate() {
        let states = [
            JointState {
                frequency1: 0.09,
                frequency2: 0.09,
                mass: 0.5,
            },
            JointState {
                frequency1: 0.0000005,
                frequency2: 0.0000005,
                mass: 0.25,
            },
            JointState {
                frequency1: 1.0,
                frequency2: 1.0,
                mass: 0.25,
            },
        ];
        let samples = [SampleConfiguration {
            acquisition: Acquisition::Poisson,
            depth: 1_000_000.0,
        }; 2];
        let (mask, certificate) = certify(&states, samples, [90_000, 90_000]).unwrap();
        assert_eq!(mask, vec![false, true, true]);
        assert!(certificate.retained_evidence_raw_log2_lower.unwrap() >= -24);
        assert!(certificate.posterior_omitted_mass_log2_upper.unwrap() < -1_000_000);
        assert_eq!(certificate.requested_mass_log2_budget, -80);
        assert_eq!(certificate.requested_mean_log2_budget, -64);
    }
    #[test]
    fn nearby_unsafe_state_with_material_mass_still_refuses_omission() {
        let states = [
            JointState {
                frequency1: 0.70365,
                frequency2: 0.70365,
                mass: 0.5,
            },
            JointState {
                frequency1: 0.7037,
                frequency2: 0.7037,
                mass: 0.5,
            },
        ];
        let samples = [SampleConfiguration {
            acquisition: Acquisition::Poisson,
            depth: 200_000.0,
        }; 2];
        assert_eq!(
            certify(&states, samples, [1, 0]).unwrap_err().reason,
            Some(CertificationReason::CumulativeOmissionBoundTooLoose)
        );
    }
}
