// SPDX-License-Identifier: GPL-3.0-or-later
//! Independent acquisition distributions. See docs/METHODS.md.
use crate::numerics::poisson_log_pmf;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AcquisitionError {
    InvalidMean,
    InvalidDispersion,
    InvalidNbpParameters,
    InvalidFrequency,
    InvalidSampleDepth,
    NonfiniteParameter(&'static str),
    DimensionMismatch,
    NumericalOverflow,
    AllocationFailure,
}
impl fmt::Display for AcquisitionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "acquisition error: {self:?}")
    }
}
impl std::error::Error for AcquisitionError {}
type Result<T> = std::result::Result<T, AcquisitionError>;
fn nonnegative(x: f64, name: &'static str, error: AcquisitionError) -> Result<()> {
    if !x.is_finite() {
        Err(AcquisitionError::NonfiniteParameter(name))
    } else if x < 0.0 {
        Err(error)
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(AcquisitionError::NumericalOverflow)
    }
}
fn positive_product(a: f64, b: f64) -> Result<f64> {
    let value = finite(a * b)?;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        Err(AcquisitionError::NumericalOverflow)
    } else {
        Ok(value)
    }
}
fn log_result(x: f64) -> Result<f64> {
    if x.is_finite() && x <= 0.0 {
        Ok(x)
    } else {
        Err(AcquisitionError::NumericalOverflow)
    }
}
fn add_logs(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let hi = a.max(b);
    hi + (-(a - b).abs()).exp().ln_1p()
}
fn softplus(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}
fn correction(x: f64) -> f64 {
    let z = x.recip();
    let q = z * z;
    z * (1.0 / 12.0 + q * (-1.0 / 360.0 + q * (1.0 / 1260.0 + q * (-1.0 / 1680.0 + q / 1188.0))))
}
fn log_factorial(k: u64) -> f64 {
    if k < 16 {
        (1..=k).map(|i| (i as f64).ln()).sum()
    } else {
        let n = k as f64;
        (n + 0.5) * n.ln() - n + 0.5 * std::f64::consts::TAU.ln() + correction(n)
    }
}
fn log_gamma_small(x: f64) -> f64 {
    let mut y = x;
    let mut shift = 0.0;
    while y < 16.0 {
        shift += y.ln();
        y += 1.0;
    }
    (y - 0.5) * y.ln() - y + 0.5 * std::f64::consts::TAU.ln() + correction(y) - shift
}
// Algebraic delta is retained even when x and y round to the same huge shape.
fn deviance(x: f64, y: f64, delta: f64) -> f64 {
    if delta.abs() <= 0.1 * x.max(y) {
        let t = (delta / x) / (1.0 + y / x);
        let mut power = t;
        let mut s = delta * t;
        for j in 1..=40 {
            power *= t * t;
            s += 2.0 * x * power / f64::from(2 * j + 1);
        }
        s
    } else {
        x * (x.ln() - y.ln()) - delta
    }
}

/// NB variance = mean + dispersion * mean^2; dispersion zero is Poisson.
#[derive(Debug, Clone, Copy)]
pub struct NegativeBinomial {
    mean: f64,
    dispersion: f64,
}
impl NegativeBinomial {
    pub fn new(mean: f64, dispersion: f64) -> Result<Self> {
        nonnegative(mean, "mean", AcquisitionError::InvalidMean)?;
        nonnegative(
            dispersion,
            "dispersion",
            AcquisitionError::InvalidDispersion,
        )?;
        Ok(Self { mean, dispersion })
    }
    pub fn log_pmf(&self, k: u64) -> Result<f64> {
        let mu = self.mean;
        let phi = self.dispersion;
        if mu == 0.0 {
            return Ok(if k == 0 { 0.0 } else { f64::NEG_INFINITY });
        }
        // Exact identity bounds the absolute Poisson-relative log correction.
        let log_bound =
            phi.ln() + 2.0 * add_logs((k as f64).ln(), mu.ln()) - std::f64::consts::LN_2;
        if phi == 0.0 || log_bound <= 1e-14_f64.ln() {
            return poisson_log_pmf(k, mu).map_err(|_| AcquisitionError::NumericalOverflow);
        }
        let log_theta = mu.ln() + phi.ln();
        let log_d = softplus(log_theta);
        let base = if log_theta < -20.0 {
            let t = log_theta.exp();
            -mu * if t == 0.0 { 1.0 } else { t.ln_1p() / t }
        } else {
            -log_d / phi
        };
        if k == 0 {
            return log_result(base);
        }
        let n = k as f64;
        let value = if k < 16 {
            let rising: f64 = (0..k).map(|i| softplus((i as f64).ln() + phi.ln())).sum();
            base + n * (mu.ln() - log_d) + rising - log_factorial(k)
        } else {
            let r = finite(phi.recip())?;
            if r == 0.0 {
                return Err(AcquisitionError::NumericalOverflow);
            }
            if r < 16.0 {
                let ratio = (r - 1.0) * n.ln() + (n + r - 0.5) * (r / n).ln_1p() - r
                    + correction(n + r)
                    - correction(n);
                ratio - log_gamma_small(r) + base - n * softplus(-log_theta)
            } else {
                let p = (-log_d).exp(); // r/(r+mu)
                let q = (-softplus(-log_theta)).exp();
                let delta_r = (mu - n) * p;
                let delta_n = (n - mu) * p;
                let target_r = r * p + n * p;
                let target_n = mu * p + n * q;
                correction(r + n)
                    - correction(r)
                    - correction(n)
                    - deviance(r, target_r, delta_r)
                    - deviance(n, target_n, delta_n)
                    - 0.5 * ((n / r).ln_1p() + std::f64::consts::TAU.ln() + n.ln())
            }
        };
        log_result(value)
    }
    pub fn pmf(&self, k: u64) -> Result<f64> {
        Ok(self.log_pmf(k)?.exp())
    }
}

/// Upstream variance law mean + beta*mean^alpha; beta and alpha are linear.
#[derive(Debug, Clone, Copy)]
pub struct VarianceLaw {
    beta: f64,
    alpha: f64,
}
impl VarianceLaw {
    pub fn new(beta: f64, alpha: f64) -> Result<Self> {
        nonnegative(beta, "beta", AcquisitionError::InvalidDispersion)?;
        if !alpha.is_finite() {
            return Err(AcquisitionError::NonfiniteParameter("alpha"));
        }
        Ok(Self { beta, alpha })
    }
    fn distribution(self, mean: f64) -> Result<NegativeBinomial> {
        if mean == 0.0 || self.beta == 0.0 {
            return NegativeBinomial::new(mean, 0.0);
        }
        let log_phi = finite(self.beta.ln() + (self.alpha - 2.0) * mean.ln())?;
        let phi = finite(log_phi.exp())?;
        if phi == 0.0 {
            return Err(AcquisitionError::NumericalOverflow);
        }
        NegativeBinomial::new(mean, phi)
    }
}
#[derive(Debug, Clone, Copy)]
pub struct NbpParameters {
    law: VarianceLaw,
    cell_depth: f64,
}
impl NbpParameters {
    pub fn new(law: VarianceLaw, cell_depth: f64) -> Result<Self> {
        nonnegative(
            cell_depth,
            "cell depth",
            AcquisitionError::InvalidNbpParameters,
        )?;
        if cell_depth == 0.0 {
            return Err(AcquisitionError::InvalidNbpParameters);
        }
        Ok(Self { law, cell_depth })
    }
    fn log_pmf(self, k: u64, f: f64, depth: f64) -> Result<f64> {
        if f == 0.0 || depth == 0.0 {
            return Ok(if k == 0 { 0.0 } else { f64::NEG_INFINITY });
        }
        let mu = positive_product(self.cell_depth, f)?;
        let c = finite(depth / self.cell_depth)?;
        if c == 0.0 {
            return Err(AcquisitionError::NumericalOverflow);
        }
        let nb = self.law.distribution(mu)?;
        let phi = nb.dispersion;
        let a = -(-c).exp_m1();
        let log_theta = mu.ln() + phi.ln();
        let log_da = if phi == 0.0 {
            0.0
        } else {
            softplus(log_theta + a.ln())
        };
        let base = if phi == 0.0 {
            -mu * a
        } else if log_theta + a.ln() < -20.0 {
            let t = (log_theta + a.ln()).exp();
            -mu * a * if t == 0.0 { 1.0 } else { t.ln_1p() / t }
        } else {
            -log_da / phi
        };
        if k == 0 {
            return log_result(base);
        }
        let len = usize::try_from(k)
            .ok()
            .and_then(|n| n.checked_add(1))
            .ok_or(AcquisitionError::DimensionMismatch)?;
        let mut stirling = storage(len, f64::NEG_INFINITY)?;
        stirling[0] = 0.0;
        for row in 1..len {
            for j in (1..=row).rev() {
                stirling[j] = add_logs(stirling[j] + (j as f64).ln(), stirling[j - 1]);
            }
            stirling[0] = f64::NEG_INFINITY;
        }
        let mut rising = 0.0;
        let mut total = f64::NEG_INFINITY;
        for (j, &s) in stirling.iter().enumerate().skip(1) {
            let log_factor = add_logs(mu.ln(), ((j - 1) as f64).ln() + log_theta) - c - log_da;
            rising += log_factor;
            total = add_logs(total, s + rising);
        }
        log_result(base + (k as f64) * c.ln() - log_factorial(k) + total)
    }
}
fn storage(n: usize, fill: f64) -> Result<Vec<f64>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| AcquisitionError::AllocationFailure)?;
    v.resize(n, fill);
    Ok(v)
}
#[derive(Debug, Clone, Copy)]
pub enum Acquisition {
    Poisson,
    NegativeBinomial(VarianceLaw),
    Nbp(NbpParameters),
}
impl Acquisition {
    fn validate(depth: f64, frequencies: &[f64]) -> Result<()> {
        nonnegative(depth, "sample depth", AcquisitionError::InvalidSampleDepth)?;
        for &f in frequencies {
            if !f.is_finite() {
                return Err(AcquisitionError::NonfiniteParameter("frequency"));
            }
            if !(0.0..=1.0).contains(&f) {
                return Err(AcquisitionError::InvalidFrequency);
            }
        }
        Ok(())
    }
    pub fn log_pmf(self, k: u64, f: f64, depth: f64) -> Result<f64> {
        Self::validate(depth, &[f])?;
        match self {
            Self::Poisson => poisson_log_pmf(k, positive_product(f, depth)?)
                .map_err(|_| AcquisitionError::NumericalOverflow),
            Self::NegativeBinomial(law) => {
                law.distribution(positive_product(f, depth)?)?.log_pmf(k)
            }
            Self::Nbp(p) => p.log_pmf(k, f, depth),
        }
    }
    pub fn pmf(self, k: u64, f: f64, depth: f64) -> Result<f64> {
        Ok(self.log_pmf(k, f, depth)?.exp())
    }
    pub fn log_vector(self, counts: &[u64], f: f64, depth: f64) -> Result<Vec<f64>> {
        Ok(self.log_matrix(counts, &[f], depth)?.values)
    }
    pub fn probability_vector(self, counts: &[u64], f: f64, depth: f64) -> Result<Vec<f64>> {
        Ok(self.probability_matrix(counts, &[f], depth)?.values)
    }
    pub fn log_matrix(self, counts: &[u64], frequencies: &[f64], depth: f64) -> Result<Matrix> {
        Self::validate(depth, frequencies)?;
        // Evaluate parameters even with no count categories.
        for &f in frequencies {
            self.log_pmf(0, f, depth)?;
        }
        let n = counts
            .len()
            .checked_mul(frequencies.len())
            .ok_or(AcquisitionError::DimensionMismatch)?;
        let mut values = storage(n, 0.0)?;
        for (i, &k) in counts.iter().enumerate() {
            for (j, &f) in frequencies.iter().enumerate() {
                values[i * frequencies.len() + j] = self.log_pmf(k, f, depth)?;
            }
        }
        Matrix::from_values(counts.len(), frequencies.len(), values)
    }
    pub fn probability_matrix(
        self,
        counts: &[u64],
        frequencies: &[f64],
        depth: f64,
    ) -> Result<Matrix> {
        let mut matrix = self.log_matrix(counts, frequencies, depth)?;
        for p in &mut matrix.values {
            *p = p.exp();
        }
        Ok(matrix)
    }
}
/// Row-major counts by frequencies. No implicit count category.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    rows: usize,
    columns: usize,
    values: Vec<f64>,
}
impl Matrix {
    pub fn from_values(rows: usize, columns: usize, values: Vec<f64>) -> Result<Self> {
        if rows.checked_mul(columns) != Some(values.len()) {
            return Err(AcquisitionError::DimensionMismatch);
        }
        if values.iter().any(|v| v.is_nan() || *v == f64::INFINITY) {
            return Err(AcquisitionError::NumericalOverflow);
        }
        Ok(Self {
            rows,
            columns,
            values,
        })
    }
    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn columns(&self) -> usize {
        self.columns
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        if row >= self.rows || column >= self.columns {
            None
        } else {
            Some(self.values[row * self.columns + column])
        }
    }
}
