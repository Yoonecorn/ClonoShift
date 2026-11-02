// SPDX-License-Identifier: GPL-3.0-or-later
//! Experimental C2-constrained Poisson fitting.
//!
//! Poisson total counts are sufficient for shared-frequency prior parameters.
//! The allocation factor is independent of those parameters; retaining only
//! (K1+K2,0) at combined depth preserves C2 and all likelihood differences.
//! This is a bounded local estimator, not biological calibration or model selection.
use crate::acquisition::Acquisition;
use crate::fitting::{self, FitConfig, FitOutcome, FitReport};
use crate::observation::{ObservedCountPair, PairTable};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AutoFitConfig {
    pub grid_points: usize,
    pub exponent_bounds: [f64; 2],
    pub log_minimum_bounds: Option<[f64; 2]>,
    pub starts: Vec<f64>,
    pub maximum_total_count: u64,
}
impl Default for AutoFitConfig {
    fn default() -> Self {
        Self {
            grid_points: 129,
            exponent_bounds: [-3.5, -0.5],
            log_minimum_bounds: None,
            starts: vec![-3.0, -2.0, -1.0],
            maximum_total_count: 1_000_000,
        }
    }
}
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct AutoCandidate {
    pub exponent: f64,
    pub minimum: f64,
    pub log_minimum: f64,
    pub c1: f64,
    pub c2: f64,
    pub mean_negative_log_likelihood: f64,
    pub total_negative_log_likelihood: f64,
    pub grid_points: usize,
}
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct SplitDiagnostic {
    pub probability_sample1: f64,
    pub informative_observations: u64,
    pub pearson_sum: f64,
    pub pearson_mean: Option<f64>,
    pub interpretation: &'static str,
}
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct StartDiagnostic {
    pub initial_exponent: f64,
    pub outcome: String,
    pub identification: String,
    pub exponent: Option<f64>,
    pub log_minimum: Option<f64>,
    pub evaluations: usize,
    pub invalid_trials: usize,
    pub constraint_slope: Option<f64>,
    pub profile_sides: Option<[f64; 2]>,
    pub refinement_c2: Option<f64>,
    pub refinement_mean_nll: Option<f64>,
}
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct GridDiagnostic {
    pub grid_points: usize,
    pub outcome: String,
    pub starts: Vec<StartDiagnostic>,
    pub incomplete_search: bool,
}
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct AutoFitReport {
    pub schema_version: u32,
    pub classification: &'static str,
    pub outcome: &'static str,
    pub reason: String,
    pub candidate: Option<AutoCandidate>,
    pub diagnostic_candidate: Option<AutoCandidate>,
    pub original_distinct_pairs: usize,
    pub total_distinct_pairs: usize,
    pub observations: u64,
    pub depths: [f64; 2],
    pub exponent_bounds: [f64; 2],
    pub log_minimum_bounds: [f64; 2],
    pub evaluations: usize,
    pub grids: Vec<GridDiagnostic>,
    pub split_pearson: Option<SplitDiagnostic>,
    pub cross_grid_exponent_difference: Option<f64>,
    pub cross_grid_log_minimum_difference: Option<f64>,
    pub parameter_uncertainty: &'static str,
    pub search_coverage: &'static str,
}
fn diagnostic(n: usize, r: &FitReport) -> GridDiagnostic {
    GridDiagnostic {
        grid_points: n,
        outcome: format!("{:?}", r.outcome),
        incomplete_search: r.incomplete_search,
        starts: r
            .starts
            .iter()
            .map(|s| StartDiagnostic {
                initial_exponent: s.initial_exponent,
                outcome: format!("{:?}", s.outcome),
                identification: format!("{:?}", s.identification),
                exponent: s.final_candidate.map(|c| c.exponent),
                log_minimum: s.final_candidate.map(|c| c.log_minimum),
                evaluations: s.evaluations,
                invalid_trials: s.invalid_trials,
                constraint_slope: s.constraint_slope,
                profile_sides: s.profile_sides,
                refinement_c2: s.refinement_c2,
                refinement_mean_nll: s.refinement_mean_nll,
            })
            .collect(),
    }
}
fn best(r: &FitReport) -> Option<fitting::Candidate> {
    r.candidate.or_else(|| {
        r.starts
            .iter()
            .filter_map(|s| s.final_candidate)
            .min_by(|a, b| {
                a.objective
                    .negative_log_likelihood
                    .total_cmp(&b.objective.negative_log_likelihood)
                    .then(a.exponent.total_cmp(&b.exponent))
                    .then(a.log_minimum.total_cmp(&b.log_minimum))
            })
    })
}
fn consistent(r: &FitReport) -> bool {
    let Some(w) = r.candidate else {
        return false;
    };
    r.outcome == FitOutcome::QualifiedFit
        && !r.incomplete_search
        && r.starts.iter().all(|s| {
            s.outcome == FitOutcome::QualifiedFit
                && s.invalid_trials == 0
                && s.final_candidate.is_some_and(|c| {
                    (c.exponent - w.exponent).abs() <= 1e-3
                        && (c.log_minimum - w.log_minimum).abs() <= 1e-3
                })
        })
}
/// Run the deterministic experimental estimator. Refusal returns the full report.
/// Numerical screens are explicit algorithmic tolerances, not calibrated tests.
pub fn fit(table: &PairTable, depths: [f64; 2], c: &AutoFitConfig) -> AutoFitReport {
    let sum_depth = depths[0] + depths[1];
    let bounds = c
        .log_minimum_bounds
        .unwrap_or([-(sum_depth.max(1.0).ln() + 8.0).max(12.0), -0.01]);
    let mut out = AutoFitReport {
        schema_version: 1, classification: "EXPERIMENTAL", outcome: "Refused",
        reason: String::new(), candidate: None, diagnostic_candidate: None,
        original_distinct_pairs: table.records().len(), total_distinct_pairs: 0,
        observations: table.observations(), depths, exponent_bounds: c.exponent_bounds,
        log_minimum_bounds: bounds, evaluations: 0, grids: vec![], split_pearson: None,
        cross_grid_exponent_difference: None, cross_grid_log_minimum_difference: None,
        parameter_uncertainty: "Unavailable; local profile diagnostics only; downstream inference is plug-in conditional",
        search_coverage: "Bounded local C2-root scans and deterministic multistart; global/root completeness not certified",
    };
    if depths
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0 || *v > (1_u64 << 53) as f64)
        || !sum_depth.is_finite()
        || table.observations() == 0
        || table.observations() > (1_u64 << 53)
        || c.grid_points < 65
        || c.grid_points > 257
        || c.maximum_total_count > 1_000_000
        || c.maximum_total_count == 0
    {
        out.reason = "Invalid or unsupported configuration: positive finite replicate depths <=2^53, observations <=2^53, 65..257 initial grid nodes and count limit <=1000000 required".into();
        return out;
    }
    let p = depths[0] / sum_depth;
    if !(0.0 < p && p < 1.0) {
        out.reason = "Replicate depth ratio is numerically unresolved".into();
        return out;
    }
    let mut counts = Vec::new();
    let mut pearson = 0.0;
    let mut information = 0_u64;
    for r in table.records() {
        let Some(t) = r.count1.checked_add(r.count2) else {
            out.reason = "Count overflow".into();
            return out;
        };
        if t > c.maximum_total_count {
            out.reason = format!(
                "Total count {t} exceeds the declared supported maximum {}",
                c.maximum_total_count
            );
            return out;
        }
        counts.push(ObservedCountPair {
            count1: t,
            count2: 0,
            multiplicity: r.multiplicity,
        });
        if t >= 2 {
            let delta = r.count1 as f64 - t as f64 * p;
            pearson += r.multiplicity as f64 * delta * delta / (t as f64 * p * (1.0 - p));
            information += r.multiplicity;
        }
    }
    if !pearson.is_finite() {
        out.reason = "Allocation diagnostic exceeds supported numerical range".into();
        return out;
    }
    out.split_pearson=Some(SplitDiagnostic { probability_sample1:p, informative_observations:information,
        pearson_sum:pearson, pearson_mean:(information>0).then(||pearson/information as f64),
        interpretation:"Descriptive conditional-binomial allocation diagnostic for totals >=2; expectation one under the shared-frequency Poisson law; no universal cutoff or biological calibration" });
    let totals = match PairTable::new(&counts) {
        Ok(t) => t,
        Err(e) => {
            out.reason = format!("Total grouping failed: {e}");
            return out;
        }
    };
    out.total_distinct_pairs = totals.records().len();
    let settings = FitConfig {
        exponent_bounds: c.exponent_bounds,
        log_minimum_bounds: bounds,
        starts: c.starts.clone(),
        grid_points: c.grid_points,
        ..FitConfig::default()
    };
    let coarse = match fitting::fit_poisson(&totals, [sum_depth, 0.0], &settings) {
        Ok(r) => r,
        Err(e) => {
            out.reason = format!("Fit configuration or observations rejected: {e:?}");
            return out;
        }
    };
    out.grids.push(diagnostic(c.grid_points, &coarse));
    // A second complete fit distinguishes parameter/grid instability from a
    // fixed-candidate refinement check. The selected fine grid is also refined
    // once internally by the inherited scout.
    let mut fine_settings = settings.clone();
    fine_settings.grid_points = 2 * c.grid_points - 1;
    let fine = match fitting::fit_poisson(&totals, [sum_depth, 0.0], &fine_settings) {
        Ok(r) => r,
        Err(e) => {
            out.reason = format!("Refined fit rejected: {e:?}");
            return out;
        }
    };
    out.grids.push(diagnostic(fine_settings.grid_points, &fine));
    out.evaluations = out
        .grids
        .iter()
        .flat_map(|g| &g.starts)
        .map(|s| s.evaluations)
        .sum();
    if let Some(w) = best(&fine) {
        if let Some(v) = fitting::evaluate_fixed(
            table,
            depths,
            Acquisition::Poisson,
            w.exponent,
            w.log_minimum,
            fine_settings.grid_points,
        ) {
            out.diagnostic_candidate = Some(AutoCandidate {
                exponent: v.exponent,
                minimum: v.log_minimum.exp(),
                log_minimum: v.log_minimum,
                c1: v.objective.c1,
                c2: v.objective.c2,
                mean_negative_log_likelihood: v.objective.mean_negative_log_likelihood,
                total_negative_log_likelihood: v.objective.negative_log_likelihood,
                grid_points: fine_settings.grid_points,
            });
        }
    }
    if let (Some(a), Some(b)) = (best(&coarse), best(&fine)) {
        out.cross_grid_exponent_difference = Some((a.exponent - b.exponent).abs());
        out.cross_grid_log_minimum_difference = Some((a.log_minimum - b.log_minimum).abs());
    }
    if !consistent(&coarse) || !consistent(&fine) {
        out.reason="No stable locally identified candidate across every starting value and both complete grid fits; see grid/start diagnostics".into();
        return out;
    }
    if out.cross_grid_exponent_difference.is_none_or(|v| v > 0.02)
        || out
            .cross_grid_log_minimum_difference
            .is_none_or(|v| v > 0.02)
    {
        out.reason =
            "Cross-grid fitted parameter difference exceeds the declared numerical screen 0.02"
                .into();
        return out;
    }
    let Some(candidate) = out.diagnostic_candidate.clone() else {
        out.reason = "Independent original-pair evaluation unavailable".into();
        return out;
    };
    if candidate.c2.abs() > 1.01e-9 {
        out.reason = "Original-pair evaluation does not satisfy corrected C2 tolerance".into();
        return out;
    }
    out.outcome = "Candidate";
    out.reason="Experimental local C2-hard Poisson fit passes declared multistart and grid screens; global optimality, parameter uncertainty and biological calibration remain unqualified".into();
    out.candidate = Some(candidate);
    out
}
