// SPDX-License-Identifier: GPL-3.0-or-later
//! Deterministic bounded Poisson fitting scout. Qualification is local and algorithmic.
//! NB/NBP APIs below are feasibility/sensitivity probes only. No parameter intervals.
use crate::acquisition::{Acquisition, NbpParameters, VarianceLaw};
use crate::numerics::FrequencyPrior;
use crate::observation::{
    ConditionedModel, Grid, ObjectiveSummary, ObservationError, PairTable, SampleConfiguration,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitOutcome {
    QualifiedFit,
    NoQualifiedFit,
    DiagnosticFail,
    SearchExhausted,
    CertifiedInfeasible,
    Unresolved,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Identification {
    LocalScreenPassed,
    Weak,
    ActiveBound,
    Unavailable,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitError {
    InvalidConfiguration,
    InvalidObservations,
}
#[derive(Debug, Clone)]
pub struct FitConfig {
    pub exponent_bounds: [f64; 2],
    pub log_minimum_bounds: [f64; 2],
    pub starts: Vec<f64>,
    pub grid_points: usize,
    pub scan_intervals: usize,
    pub root_iterations: usize,
    pub search_iterations: usize,
    pub equality_tolerance: f64,
    pub root_tolerance: f64,
    pub step_tolerance: f64,
    pub refinement_tolerance: f64,
    pub identification_tolerance: f64,
}
impl Default for FitConfig {
    fn default() -> Self {
        Self {
            exponent_bounds: [-3.0, 1.0],
            log_minimum_bounds: [-12.0, -0.01],
            starts: vec![-2.0, -1.0, 0.0],
            grid_points: 65,
            scan_intervals: 32,
            root_iterations: 80,
            search_iterations: 100,
            equality_tolerance: 1e-9,
            root_tolerance: 1e-11,
            step_tolerance: 1e-4,
            refinement_tolerance: 5e-4,
            identification_tolerance: 1e-7,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Candidate {
    pub exponent: f64,
    pub log_minimum: f64,
    pub objective: ObjectiveSummary,
}
#[derive(Debug, Clone, PartialEq)]
pub struct StartRecord {
    pub id: usize,
    pub initial_exponent: f64,
    pub final_candidate: Option<Candidate>,
    pub outcome: FitOutcome,
    pub identification: Identification,
    pub iterations: usize,
    pub evaluations: usize,
    pub invalid_trials: usize,
    pub constraint_slope: Option<f64>,
    pub profile_change: Option<f64>,
    pub profile_sides: Option<[f64; 2]>,
    pub evaluation_errors: std::collections::BTreeMap<String, usize>,
    pub refinement_c2: Option<f64>,
    pub refinement_mean_nll: Option<f64>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct FitReport {
    pub outcome: FitOutcome,
    pub candidate: Option<Candidate>,
    pub starts: Vec<StartRecord>,
    pub incomplete_search: bool,
    pub reason: &'static str,
}
fn valid(c: &FitConfig, depths: [f64; 2], table: &PairTable) -> Result<(), FitError> {
    let ordered = |v: [f64; 2]| v[0].is_finite() && v[1].is_finite() && v[0] < v[1];
    if !ordered(c.exponent_bounds)
        || !ordered(c.log_minimum_bounds)
        || c.log_minimum_bounds[1] >= 0.0
        || c.log_minimum_bounds[0].exp() == 0.0
        || c.log_minimum_bounds[1].exp() >= 1.0
        || !(c.log_minimum_bounds[1] - c.log_minimum_bounds[0]).is_finite()
        || !(c.exponent_bounds[1] - c.exponent_bounds[0]).is_finite()
        || c.starts.is_empty()
        || c.starts.len() > 64
        || c.starts
            .iter()
            .any(|a| !a.is_finite() || *a < c.exponent_bounds[0] || *a > c.exponent_bounds[1])
        || c.grid_points < 3
        || c.grid_points > 1025
        || c.scan_intervals < 2
        || c.scan_intervals > 4096
        || c.root_iterations > 256
        || c.search_iterations > 10000
        || [
            c.equality_tolerance,
            c.root_tolerance,
            c.step_tolerance,
            c.refinement_tolerance,
            c.identification_tolerance,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || depths.iter().any(|d| !d.is_finite() || *d < 0.0)
        || depths[0] + depths[1] <= 0.0
        || !(depths[0] + depths[1]).is_finite()
    {
        return Err(FitError::InvalidConfiguration);
    }
    // Admit the complete caller box only when the accepted prior can represent
    // its support-endpoint log densities. Log density is affine in ln(f), so
    // these corners check its extrema without an arbitrary exponent cap.
    // This is a necessary prior-arithmetic check, not objective qualification.
    for exponent in c.exponent_bounds {
        for log_minimum in c.log_minimum_bounds {
            let minimum = log_minimum.exp();
            let prior = FrequencyPrior::new(exponent, minimum)
                .map_err(|_| FitError::InvalidConfiguration)?;
            for frequency in [minimum, 1.0] {
                if !prior.log_density(frequency).is_ok_and(|v| v.is_finite()) {
                    return Err(FitError::InvalidConfiguration);
                }
            }
        }
    }
    if table.observations() == 0
        || table
            .records()
            .iter()
            .any(|r| (depths[0] == 0.0 && r.count1 > 0) || (depths[1] == 0.0 && r.count2 > 0))
    {
        return Err(FitError::InvalidObservations);
    }
    Ok(())
}
pub fn evaluate_fixed(
    table: &PairTable,
    depths: [f64; 2],
    acquisition: Acquisition,
    a: f64,
    x: f64,
    points: usize,
) -> Option<Candidate> {
    evaluate_result(table, depths, acquisition, a, x, points).ok()
}
fn evaluate_result(
    table: &PairTable,
    depths: [f64; 2],
    acquisition: Acquisition,
    a: f64,
    x: f64,
    points: usize,
) -> Result<Candidate, ObservationError> {
    let grid = Grid::log_uniform(x.exp(), points)?;
    let model = ConditionedModel::new(
        grid,
        a,
        SampleConfiguration {
            acquisition,
            depth: depths[0],
        },
        SampleConfiguration {
            acquisition,
            depth: depths[1],
        },
    )?;
    Ok(Candidate {
        exponent: a,
        log_minimum: x,
        objective: model.evaluate(table)?,
    })
}

fn better(a: Candidate, b: Candidate) -> bool {
    let d = a.objective.negative_log_likelihood - b.objective.negative_log_likelihood;
    if d.abs() > 1e-10 {
        d < 0.0
    } else {
        a.exponent
            .total_cmp(&b.exponent)
            .then(a.log_minimum.total_cmp(&b.log_minimum))
            .is_lt()
    }
}
struct Search<'a> {
    table: &'a PairTable,
    depths: [f64; 2],
    acquisition: Acquisition,
    c: &'a FitConfig,
    evaluations: usize,
    invalid: usize,
    errors: std::collections::BTreeMap<String, usize>,
    exhausted: bool,
}
impl Search<'_> {
    fn eval(&mut self, a: f64, x: f64) -> Option<Candidate> {
        self.evaluations += 1;
        let v = evaluate_result(
            self.table,
            self.depths,
            self.acquisition,
            a,
            x,
            self.c.grid_points,
        );
        match v {
            Ok(v) => Some(v),
            Err(e) => {
                self.invalid += 1;
                *self.errors.entry(format!("{e:?}")).or_default() += 1;
                None
            }
        }
    }
    fn root(&mut self, mut left: Candidate, mut right: Candidate) -> Option<Candidate> {
        for _ in 0..self.c.root_iterations {
            if left.objective.c2.abs() <= self.c.equality_tolerance {
                return Some(left);
            }
            if right.objective.c2.abs() <= self.c.equality_tolerance {
                return Some(right);
            }
            let middle = self.eval(left.exponent, (left.log_minimum + right.log_minimum) / 2.0)?;
            if middle.objective.c2.abs() <= self.c.equality_tolerance {
                return Some(middle);
            }
            if right.log_minimum - left.log_minimum <= self.c.root_tolerance {
                // Width convergence alone never establishes equality feasibility.
                self.exhausted = true;
                return None;
            }
            if left.objective.c2.signum() != middle.objective.c2.signum() {
                right = middle;
            } else {
                left = middle;
            }
        }
        self.exhausted = true;
        None
    }
    fn feasible(&mut self, a: f64) -> Option<Candidate> {
        let mut previous = None;
        let mut best = None;
        for i in 0..=self.c.scan_intervals {
            let x = self.c.log_minimum_bounds[0]
                + (self.c.log_minimum_bounds[1] - self.c.log_minimum_bounds[0]) * i as f64
                    / self.c.scan_intervals as f64;
            let current = self.eval(a, x);
            if let Some(p) = current {
                let root = if p.objective.c2.abs() <= self.c.equality_tolerance {
                    Some(p)
                } else if let Some(q) = previous {
                    let q: Candidate = q;
                    if q.objective.c2.signum() != p.objective.c2.signum() {
                        self.root(q, p)
                    } else {
                        None
                    }
                } else {
                    None
                };
                if let Some(r) = root {
                    if best.is_none_or(|b| better(r, b)) {
                        best = Some(r);
                    }
                }
            }
            previous = current; // Never bracket across invalid evaluations.
        }
        best
    }
}
/// Deterministic multi-start reduced search; all intervals are Unavailable.
/// No C1 threshold exists. DiagnosticFail is reserved for later calibrated diagnostics.
pub fn fit_poisson(
    table: &PairTable,
    depths: [f64; 2],
    c: &FitConfig,
) -> Result<FitReport, FitError> {
    valid(c, depths, table)?;
    if (table.observations() as f64) * c.log_minimum_bounds[0].exp() > 1.0 {
        return Ok(FitReport {
            outcome: FitOutcome::CertifiedInfeasible,
            candidate: None,
            starts: vec![],
            incomplete_search: false,
            reason: "C2 observed mass >= Nobs*minimum_frequency > 1 throughout domain",
        });
    }
    let mut records = Vec::new();
    let mut winner = None;
    let mut incomplete = false;
    for (id, &initial) in c.starts.iter().enumerate() {
        let mut s = Search {
            table,
            depths,
            acquisition: Acquisition::Poisson,
            c,
            evaluations: 0,
            invalid: 0,
            errors: std::collections::BTreeMap::new(),
            exhausted: false,
        };
        let mut a = initial;
        let mut best = s.feasible(a);
        let mut step = (c.exponent_bounds[1] - c.exponent_bounds[0]) / 4.0;
        let mut converged = false;
        let mut iterations = 0;
        for iteration in 0..c.search_iterations {
            iterations = iteration + 1;
            let mut next = best;
            for trial in [
                (a - step).max(c.exponent_bounds[0]),
                (a + step).min(c.exponent_bounds[1]),
            ] {
                if let Some(v) = s.feasible(trial) {
                    if next.is_none_or(|b| better(v, b)) {
                        next = Some(v);
                    }
                }
            }
            if next != best {
                best = next;
                a = best.expect("candidate exists").exponent;
            } else {
                step /= 2.0;
            }
            if step <= c.step_tolerance {
                converged = true;
                break;
            }
        }
        let mut identification = Identification::Unavailable;
        let mut profile_sides = None;
        let (mut slope, mut change, mut refine_c2, mut refine_nll) = (None, None, None, None);
        let mut outcome = if s.invalid == s.evaluations {
            FitOutcome::Unresolved
        } else if !converged || s.exhausted {
            FitOutcome::SearchExhausted
        } else if s.invalid > 0 {
            FitOutcome::Unresolved
        } else {
            FitOutcome::NoQualifiedFit
        };
        if let Some(b) = best {
            if converged && !s.exhausted {
                let h = 1e-4;
                if let (Some(l), Some(r)) = (
                    s.eval(b.exponent, b.log_minimum - h),
                    s.eval(b.exponent, b.log_minimum + h),
                ) {
                    slope = Some((r.objective.c2 - l.objective.c2) / (2.0 * h));
                }
                let mut differences = Vec::new();
                for t in [
                    (b.exponent - 0.05).max(c.exponent_bounds[0]),
                    (b.exponent + 0.05).min(c.exponent_bounds[1]),
                ] {
                    if let Some(p) = s.feasible(t) {
                        differences.push(
                            p.objective.mean_negative_log_likelihood
                                - b.objective.mean_negative_log_likelihood,
                        );
                    }
                }
                if differences.len() == 2 {
                    profile_sides = Some([differences[0], differences[1]]);
                    change = Some(differences.into_iter().fold(f64::INFINITY, f64::min));
                }
                if let Some(f) = evaluate_fixed(
                    table,
                    depths,
                    Acquisition::Poisson,
                    b.exponent,
                    b.log_minimum,
                    2 * c.grid_points - 1,
                ) {
                    refine_c2 = Some(f.objective.c2);
                    refine_nll = Some(
                        f.objective.mean_negative_log_likelihood
                            - b.objective.mean_negative_log_likelihood,
                    );
                }
                let bound = (b.exponent - c.exponent_bounds[0]).abs() <= 2.0 * c.step_tolerance
                    || (b.exponent - c.exponent_bounds[1]).abs() <= 2.0 * c.step_tolerance
                    || (b.log_minimum - c.log_minimum_bounds[0]).abs() <= 2.0 * c.root_tolerance
                    || (b.log_minimum - c.log_minimum_bounds[1]).abs() <= 2.0 * c.root_tolerance;
                identification = if bound {
                    Identification::ActiveBound
                } else if slope.is_none_or(|v| v.abs() < 1e-6)
                    || change.is_none_or(|v| v < c.identification_tolerance)
                {
                    Identification::Weak
                } else {
                    Identification::LocalScreenPassed
                };
                outcome = if s.exhausted {
                    FitOutcome::SearchExhausted
                } else if b.objective.c2.abs() > c.equality_tolerance {
                    FitOutcome::NoQualifiedFit
                } else if identification != Identification::LocalScreenPassed
                    || refine_c2.is_none_or(|v| v.abs() > c.refinement_tolerance)
                    || refine_nll.is_none_or(|v| v.abs() > c.refinement_tolerance)
                {
                    FitOutcome::Unresolved
                } else {
                    FitOutcome::QualifiedFit
                };
            }
        }
        incomplete |= outcome == FitOutcome::SearchExhausted;
        if outcome == FitOutcome::QualifiedFit {
            let b = best.unwrap();
            if winner.is_none_or(|w| better(b, w)) {
                winner = Some(b);
            }
        }
        records.push(StartRecord {
            id,
            initial_exponent: initial,
            final_candidate: best,
            outcome,
            identification,
            iterations,
            evaluations: s.evaluations,
            invalid_trials: s.invalid,
            constraint_slope: slope,
            profile_change: change,
            profile_sides,
            evaluation_errors: s.errors,
            refinement_c2: refine_c2,
            refinement_mean_nll: refine_nll,
        });
    }
    let outcome = if incomplete {
        FitOutcome::SearchExhausted
    } else if winner.is_some() {
        FitOutcome::QualifiedFit
    } else if records.iter().any(|r| r.outcome == FitOutcome::Unresolved) {
        FitOutcome::Unresolved
    } else {
        FitOutcome::NoQualifiedFit
    };
    Ok(FitReport {outcome,candidate:winner,starts:records,incomplete_search:incomplete,reason:"local scout only; root scan/global optimum and parameter intervals unqualified; C1 diagnostic uncalibrated"})
}
#[derive(Debug, Clone, PartialEq)]
pub struct Comparator {
    pub outcome: FitOutcome,
    pub candidate: Option<Candidate>,
    pub evaluations: usize,
    pub incomplete_search: bool,
}
/// Both-hard feasibility on the minimum-NLL discovered C2 branch; branch coverage incomplete.
pub fn both_hard_comparator(
    table: &PairTable,
    depths: [f64; 2],
    c: &FitConfig,
) -> Result<Comparator, FitError> {
    valid(c, depths, table)?;
    if depths[0] + depths[1] <= table.observations() as f64 {
        return Ok(Comparator {
            outcome: FitOutcome::CertifiedInfeasible,
            candidate: None,
            evaluations: 0,
            incomplete_search: false,
        });
    }
    let mut s = Search {
        table,
        depths,
        acquisition: Acquisition::Poisson,
        c,
        evaluations: 0,
        invalid: 0,
        errors: std::collections::BTreeMap::new(),
        exhausted: false,
    };
    let mut prev: Option<Candidate> = None;
    let mut best = None;
    for i in 0..=c.scan_intervals {
        let a = c.exponent_bounds[0]
            + (c.exponent_bounds[1] - c.exponent_bounds[0]) * i as f64 / c.scan_intervals as f64;
        let current = s.feasible(a);
        if let Some(v) = current {
            if v.objective.c1.abs() <= c.equality_tolerance {
                best = Some(v);
            }
            if let Some(mut l) = prev {
                let mut r = v;
                if l.objective.c1.signum() != r.objective.c1.signum() {
                    let mut found = false;
                    for _ in 0..c.root_iterations {
                        if let Some(m) = s.feasible((l.exponent + r.exponent) / 2.0) {
                            if m.objective.c1.abs() <= c.equality_tolerance {
                                best = Some(m);
                                found = true;
                                break;
                            }
                            if l.objective.c1.signum() != m.objective.c1.signum() {
                                r = m;
                            } else {
                                l = m;
                            }
                        } else {
                            break;
                        }
                    }
                    if !found {
                        s.exhausted = true;
                    }
                }
            }
        }
        prev = current;
    }
    Ok(Comparator {
        outcome: if s.exhausted {
            FitOutcome::SearchExhausted
        } else if best.is_some() || s.invalid > 0 {
            FitOutcome::Unresolved
        } else {
            FitOutcome::NoQualifiedFit
        },
        candidate: best,
        evaluations: s.evaluations,
        incomplete_search: true,
    })
}
/// Bounded unconstrained coordinate pattern search, solely a reference comparator.
pub fn unconstrained_reference(
    table: &PairTable,
    depths: [f64; 2],
    c: &FitConfig,
) -> Result<Comparator, FitError> {
    valid(c, depths, table)?;
    let mut best = None;
    let mut evaluations = 0;
    let mut exhausted = false;
    for &a in &c.starts {
        let mut local = evaluate_fixed(
            table,
            depths,
            Acquisition::Poisson,
            a,
            (c.log_minimum_bounds[0] + c.log_minimum_bounds[1]) / 2.0,
            c.grid_points,
        );
        evaluations += 1;
        let mut steps = [
            (c.exponent_bounds[1] - c.exponent_bounds[0]) / 4.0,
            (c.log_minimum_bounds[1] - c.log_minimum_bounds[0]) / 4.0,
        ];
        let mut done = false;
        for _ in 0..c.search_iterations {
            let Some(b) = local else {
                break;
            };
            let mut next = local;
            for (a, x) in [
                (b.exponent - steps[0], b.log_minimum),
                (b.exponent + steps[0], b.log_minimum),
                (b.exponent, b.log_minimum - steps[1]),
                (b.exponent, b.log_minimum + steps[1]),
            ] {
                evaluations += 1;
                if let Some(v) = evaluate_fixed(
                    table,
                    depths,
                    Acquisition::Poisson,
                    a.clamp(c.exponent_bounds[0], c.exponent_bounds[1]),
                    x.clamp(c.log_minimum_bounds[0], c.log_minimum_bounds[1]),
                    c.grid_points,
                ) {
                    if next.is_none_or(|w| better(v, w)) {
                        next = Some(v);
                    }
                }
            }
            if next == local {
                steps = steps.map(|v| v / 2.0);
            }
            local = next;
            if steps[0].max(steps[1]) <= c.step_tolerance {
                done = true;
                break;
            }
        }
        exhausted |= !done;
        if let Some(b) = local {
            if best.is_none_or(|w| better(b, w)) {
                best = Some(b);
            }
        }
    }
    Ok(Comparator {
        outcome: if exhausted {
            FitOutcome::SearchExhausted
        } else {
            FitOutcome::Unresolved
        },
        candidate: best,
        evaluations,
        incomplete_search: exhausted,
    })
}
/// Fixed-acquisition C2 root feasibility; NB/NBP are never promoted to QualifiedFit.
pub fn feasibility_profile(
    table: &PairTable,
    depths: [f64; 2],
    model: Acquisition,
    c: &FitConfig,
) -> Result<Comparator, FitError> {
    valid(c, depths, table)?;
    let mut s = Search {
        table,
        depths,
        acquisition: model,
        c,
        evaluations: 0,
        invalid: 0,
        errors: std::collections::BTreeMap::new(),
        exhausted: false,
    };
    let mut best = None;
    for &a in &c.starts {
        if let Some(b) = s.feasible(a) {
            if best.is_none_or(|w| better(b, w)) {
                best = Some(b);
            }
        }
    }
    Ok(Comparator {
        outcome: if s.exhausted {
            FitOutcome::SearchExhausted
        } else if best.is_some() || s.invalid > 0 {
            FitOutcome::Unresolved
        } else {
            FitOutcome::NoQualifiedFit
        },
        candidate: best,
        evaluations: s.evaluations,
        incomplete_search: true,
    })
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sensitivity {
    pub alpha_unidentified: bool,
    pub m_sensitive: bool,
    pub beta_m_cosine: Option<f64>,
    pub weak_beta_m: bool,
    pub m_norm: f64,
}
/// Local finite-difference prediction screen on caller frequencies/count categories.
/// Cosine magnitude >0.995 or insensitive M flags weak separation. No interval inferred.
pub fn nbp_sensitivity(
    beta: f64,
    alpha: f64,
    m: f64,
    depth: f64,
    frequencies: &[f64],
    counts: &[u64],
) -> Result<Sensitivity, FitError> {
    let make = |b, mm| {
        NbpParameters::new(VarianceLaw::new(b, alpha).ok()?, mm)
            .ok()
            .map(Acquisition::Nbp)
    };
    if frequencies.is_empty() || counts.is_empty() || !m.is_finite() || m <= 0.0 {
        return Err(FitError::InvalidConfiguration);
    }
    let base = make(beta, m).ok_or(FitError::InvalidConfiguration)?;
    let h: f64 = 1e-3;
    let db = 1e-4 * (1.0 + beta);
    let bm = make((beta - db).max(0.0), m).ok_or(FitError::InvalidConfiguration)?;
    let bp = make(beta + db, m).ok_or(FitError::InvalidConfiguration)?;
    let mm = make(beta, m * (-h).exp()).ok_or(FitError::InvalidConfiguration)?;
    let mp = make(beta, m * h.exp()).ok_or(FitError::InvalidConfiguration)?;
    let (mut bb, mut vv, mut bv) = (0.0, 0.0, 0.0);
    for &f in frequencies {
        for &k in counts {
            let p = base
                .log_pmf(k, f, depth)
                .map_err(|_| FitError::InvalidConfiguration)?;
            let b = (bp
                .log_pmf(k, f, depth)
                .map_err(|_| FitError::InvalidConfiguration)?
                - bm.log_pmf(k, f, depth)
                    .map_err(|_| FitError::InvalidConfiguration)?)
                / (beta + db - (beta - db).max(0.0));
            let v = (mp
                .log_pmf(k, f, depth)
                .map_err(|_| FitError::InvalidConfiguration)?
                - mm.log_pmf(k, f, depth)
                    .map_err(|_| FitError::InvalidConfiguration)?)
                / (2.0 * h);
            if ![p, b, v].iter().all(|v| v.is_finite()) {
                return Err(FitError::InvalidConfiguration);
            }
            // Probability-weighted prediction sensitivity, not a fitted Fisher matrix.
            bb += p.exp() * b * b;
            vv += p.exp() * v * v;
            bv += p.exp() * b * v;
        }
    }
    if ![bb, vv, bv].iter().all(|v| v.is_finite()) {
        return Err(FitError::InvalidConfiguration);
    }
    let cosine = if bb > 0.0 && vv > 0.0 {
        let cos = (bv / bb.sqrt()) / vv.sqrt();
        if !cos.is_finite() || cos.abs() > 1.0 + 64.0 * f64::EPSILON {
            return Err(FitError::InvalidConfiguration);
        }
        Some(cos)
    } else {
        None
    };
    let sensitive = vv.sqrt() > 1e-6;
    Ok(Sensitivity {
        alpha_unidentified: beta == 0.0,
        m_sensitive: sensitive,
        beta_m_cosine: cosine,
        weak_beta_m: !sensitive || cosine.is_none_or(|v| v.abs() > 0.995),
        m_norm: vv.sqrt(),
    })
}
