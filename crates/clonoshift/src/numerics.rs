// SPDX-License-Identifier: GPL-3.0-or-later
//! Independent f64 Poisson and normalized truncated power-law kernels.
//! No conditioned likelihood, inference, or fitting is implemented here.

use std::fmt;

/// Invalid scientific domains are distinguished from finite arithmetic range limits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumericalError {
    Domain(&'static str),
    NumericalRange(&'static str),
}
impl fmt::Display for NumericalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Domain(message) => write!(f, "domain error: {message}"),
            Self::NumericalRange(message) => write!(f, "numerical range error: {message}"),
        }
    }
}
impl std::error::Error for NumericalError {}
type Result<T> = std::result::Result<T, NumericalError>;

fn finite(value: f64, message: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NumericalError::NumericalRange(message))
    }
}

/// Log P(k|lambda) = k ln(lambda) - lambda - ln(k!).
///
/// k is unsigned; lambda must be finite and nonnegative. At lambda=0,
/// k=0 returns 0, otherwise negative infinity (mathematically zero mass).
/// Stable deviance/Stirling evaluation avoids factorial overflow and mode cancellation.
/// Counts above 2^53 undergo the documented u64-to-f64 rounding.
pub fn poisson_log_pmf(k: u64, lambda: f64) -> Result<f64> {
    if !lambda.is_finite() || lambda < 0.0 {
        return Err(NumericalError::Domain(
            "Poisson mean must be finite and nonnegative",
        ));
    }
    if k == 0 {
        return Ok(-lambda);
    }
    if lambda == 0.0 {
        return Ok(f64::NEG_INFINITY);
    }
    let n = k as f64;
    let value = if k < 16 {
        let log_factorial: f64 = (1..=k).map(|j| (j as f64).ln()).sum();
        n * lambda.ln() - lambda - log_factorial
    } else {
        let inverse = n.recip();
        let squared = inverse * inverse;
        let stirling = inverse
            * (1.0 / 12.0
                + squared
                    * (-1.0 / 360.0
                        + squared * (1.0 / 1260.0 + squared * (-1.0 / 1680.0 + squared / 1188.0))));
        let deviance = if lambda >= n / 2.0 && lambda <= n * 2.0 {
            let delta = n - lambda;
            let t = delta / (n + lambda);
            let mut power = t;
            let mut correction = 0.0;
            for j in 1..=40 {
                power *= t * t;
                correction += power / f64::from(2 * j + 1);
            }
            delta * t + 2.0 * n * correction
        } else if n > lambda {
            n * (n.ln() - lambda.ln()) + (lambda - n)
        } else {
            (lambda - n) - n * (lambda.ln() - n.ln())
        };
        -deviance - 0.5 * (std::f64::consts::TAU.ln() + n.ln()) - stirling
    };
    let value = finite(value, "positive-mean Poisson log-PMF is not representable")?;
    if value > 0.0 {
        return Err(NumericalError::NumericalRange(
            "Poisson log-PMF exceeded zero",
        ));
    }
    Ok(value)
}

/// Direct Poisson mass. Underflow may yield zero; the finite log-PMF is authoritative.
pub fn poisson_pmf(k: u64, lambda: f64) -> Result<f64> {
    let value = poisson_log_pmf(k, lambda)?.exp();
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(NumericalError::NumericalRange(
            "Poisson PMF is outside [0,1]",
        ))
    }
}

/// Continuous density rho(f)=f^exponent / integral_minimum^1 f^exponent df.
///
/// Finite exponent and 0<minimum<1 are required. Support [minimum,1] is closed.
/// Normalization is analytical, intentionally distinct from NoisET's grid trapezoid.
#[derive(Debug, Clone, Copy)]
pub struct FrequencyPrior {
    exponent: f64,
    minimum: f64,
    log_width: f64,
    log_scale: f64,
}
impl FrequencyPrior {
    pub fn new(exponent: f64, minimum: f64) -> Result<Self> {
        if !exponent.is_finite() || !minimum.is_finite() || minimum <= 0.0 || minimum >= 1.0 {
            return Err(NumericalError::Domain(
                "prior requires finite exponent and finite 0<minimum<1",
            ));
        }
        let log_width = -minimum.ln();
        let b = exponent + 1.0;
        let log_scale = if b == 0.0 {
            -log_width.ln()
        } else {
            // If this product overflows, exp(-product) is necessarily below f64
            // representation; exp_m1(-infinity)=-1 gives the analytic correction.
            let correction = (-(-b.abs() * log_width).exp_m1()).ln();
            b.abs().ln() - correction
        };
        finite(log_scale, "prior normalization scale is not representable")?;
        Ok(Self {
            exponent,
            minimum,
            log_width,
            log_scale,
        })
    }

    /// Log density on the closed support; finite outside-support f gives -infinity.
    /// NaN/infinite coordinates are Domain errors. Supported log overflow is a
    /// NumericalRange error, never a fictitious zero density or NaN sentinel.
    pub fn log_density(&self, f: f64) -> Result<f64> {
        if !f.is_finite() {
            return Err(NumericalError::Domain("frequency must be finite"));
        }
        if f < self.minimum || f > 1.0 {
            return Ok(f64::NEG_INFINITY);
        }
        let b = self.exponent + 1.0;
        let value = if b == 0.0 {
            -f.ln() + self.log_scale
        } else if b > 0.0 {
            self.exponent * f.ln() + self.log_scale
        } else {
            let distance = if f <= self.minimum * 2.0 {
                ((f - self.minimum) / self.minimum).ln_1p()
            } else {
                f.ln() + self.log_width
            };
            self.exponent * distance + self.log_width + self.log_scale
        };
        finite(value, "supported prior log density is not representable")
    }

    /// Direct density may underflow to zero. Overflow is NumericalRange; density
    /// is not a probability mass and is allowed to exceed one.
    pub fn density(&self, f: f64) -> Result<f64> {
        finite(
            self.log_density(f)?.exp(),
            "prior density overflow; use log_density",
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_poisson_boundaries_and_domains() {
        assert_eq!(poisson_log_pmf(0, 0.0).unwrap(), 0.0);
        assert_eq!(poisson_pmf(0, 0.0).unwrap(), 1.0);
        assert_eq!(poisson_log_pmf(1, 0.0).unwrap(), f64::NEG_INFINITY);
        assert_eq!(poisson_pmf(1, 0.0).unwrap(), 0.0);
        for invalid in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                poisson_log_pmf(1, invalid),
                Err(NumericalError::Domain(_))
            ));
            assert!(matches!(
                poisson_pmf(1, invalid),
                Err(NumericalError::Domain(_))
            ));
        }
        assert_eq!(poisson_log_pmf(0, 1.0).unwrap(), -1.0);
        assert_eq!(poisson_log_pmf(1, 1.0).unwrap(), -1.0);
    }
    #[test]
    fn prior_support_domains_and_range_errors() {
        let p = FrequencyPrior::new(-1.0, 0.25).unwrap();
        for f in [0.0, -1.0, 0.125, 2.0] {
            assert_eq!(p.log_density(f).unwrap(), f64::NEG_INFINITY);
            assert_eq!(p.density(f).unwrap(), 0.0);
        }
        for f in [0.25, 0.5, 1.0] {
            assert!(p.log_density(f).unwrap().is_finite());
        }
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(p.log_density(x), Err(NumericalError::Domain(_))));
            assert!(matches!(p.density(x), Err(NumericalError::Domain(_))));
            assert!(FrequencyPrior::new(x, 0.25).is_err());
        }
        for cutoff in [0.0, -0.1, 1.0, 2.0, f64::NAN, f64::INFINITY] {
            assert!(FrequencyPrior::new(0.0, cutoff).is_err());
        }
        let steep = FrequencyPrior::new(-f64::MAX, f64::from_bits(1)).unwrap();
        assert!(steep.log_density(f64::from_bits(1)).unwrap().is_finite());
        assert!(matches!(
            steep.density(f64::from_bits(1)),
            Err(NumericalError::NumericalRange(_))
        ));
        assert!(matches!(
            steep.log_density(1.0),
            Err(NumericalError::NumericalRange(_))
        ));
        let positive = FrequencyPrior::new(f64::MAX, 0.5).unwrap();
        assert!(positive.log_density(1.0).unwrap().is_finite());
    }
}
