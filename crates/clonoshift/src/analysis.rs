// SPDX-License-Identifier: GPL-3.0-or-later
//! Fixed supplied-parameter integration. Existing scientific components remain
//! authoritative; no fitting, model selection or biological decision rule.
use crate::acquisition::Acquisition;
use crate::observation::{ConditionedModel, ObservedCountPair, PairTable, SampleConfiguration};
use crate::response::{ResponseError, ResponseModel};
use crate::response_prior::{build_symmetric_prior, Baseline, ResponsePriorConfig};
use serde::ser::{SerializeMap, SerializeSeq};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::io::{self, BufWriter, Write};
use std::{collections::BTreeMap, fmt, fs, path::Path};

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum State {
    InvalidInput,
    InvalidConfiguration,
    UnsupportedOperation,
    NumericalFailure,
    UnresolvedResult,
    LegacyReferenceFailure,
}
#[derive(Debug)]
pub struct AnalysisError {
    pub state: State,
    pub message: String,
    pub certification_reason: Option<crate::certification::CertificationReason>,
}
impl fmt::Display for AnalysisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.state, self.message)
    }
}
impl std::error::Error for AnalysisError {}
type Result<T> = std::result::Result<T, AnalysisError>;
fn err(state: State, message: impl ToString) -> AnalysisError {
    AnalysisError {
        state,
        message: message.to_string(),
        certification_reason: None,
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Noise {
    pub model: String,
    pub parameter_order: [String; 2],
    /// Exponent then minimum frequency (linear, never log10).
    pub parameters: [f64; 2],
    pub source: String,
    pub source_reference: String,
    #[serde(default)]
    pub legacy_model_file: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Prior {
    pub response_probability: f64,
    pub scale: f64,
    pub points: usize,
    pub stride: usize,
    pub maximum_shift: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub noise: Noise,
    pub response_prior: Prior,
    pub effect_delta: f64,
    #[serde(default)]
    pub second_baseline: Option<[f64; 2]>,
    #[serde(default = "fixed_operation")]
    pub operation: String,
}
fn fixed_operation() -> String {
    "fixed_parameter".into()
}
#[derive(Debug, Clone)]
pub struct Analysis {
    pub metadata: Value,
    pub data: Value,
    pub comparison: Option<Value>,
}

// Internal file-workflow records use fixed storage. The public Value-based
// Analysis API remains an explicit conversion, with no separate inference path.
#[derive(Serialize)]
struct Clonotype {
    #[serde(rename = "AACDR3")]
    aa: String,
    counts: [u64; 2],
    credible_interval_95: [f64; 2],
    log_evidence: f64,
    negative_shift_probability: f64,
    no_qualified_change_probability: f64,
    #[serde(rename = "ntCDR3")]
    nt: String,
    null_pair: NullPair,
    positive_shift_probability: f64,
    posterior_mass: f64,
    posterior_response_mean: f64,
    posterior_underflows: usize,
    probability_within_delta: f64,
    scientific_directional_call: Option<()>,
    state: &'static str,
}
#[derive(Serialize)]
struct NullPair {
    log_conditioned: f64,
    log_posterior_mean_frequency: f64,
    log_unconditioned: f64,
}
struct Prepared {
    metadata: Value,
    data: Value,
    records: Vec<Clonotype>,
    comparison: Option<PreparedComparison>,
}
struct PreparedComparison {
    header: Value,
    summaries: Vec<crate::legacy::Posterior>,
    indices: Vec<usize>,
    certified: bool,
}
// Headers retain their sorted Value keys, including the clonotypes placeholder.
// Replace only its value during serialization; no per-clone Value is built.
struct Clonotypes<'a, T> {
    header: &'a Value,
    records: T,
}
impl<T: Serialize> Serialize for Clonotypes<'_, T> {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        let header = self.header.as_object().expect("internal object header");
        let mut map = serializer.serialize_map(Some(header.len()))?;
        for (key, value) in header {
            if key == "clonotypes" {
                map.serialize_entry(key, &self.records)?;
            } else {
                map.serialize_entry(key, value)?;
            }
        }
        map.end()
    }
}
#[derive(Serialize)]
struct LegacySummary<'a> {
    lower: Option<f64>,
    mean: f64,
    median: Option<f64>,
    mode: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    posterior: Option<&'a [f64]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    shared_posterior_index: Option<usize>,
    tail: f64,
    upper: Option<f64>,
}
#[derive(Serialize)]
struct Signs {
    corrected_mean_sign: i32,
    legacy_median_sign: Option<i32>,
    signs_agree: Option<bool>,
}
#[derive(Serialize)]
struct ComparisonRecord<'a> {
    corrected_response: &'a Clonotype,
    counts: [u64; 2],
    defined_differences: &'static [&'static str],
    descriptive_signs: Signs,
    direct_scientific_comparison: &'static str,
    directional_interpretation_agrees: Option<()>,
    legacy_response: LegacySummary<'a>,
    #[serde(rename = "ntCDR3")]
    nt: &'a str,
    reason: &'static str,
}
struct ComparisonRows<'a> {
    records: &'a [Clonotype],
    comparison: &'a PreparedComparison,
}
impl Serialize for ComparisonRows<'_> {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        let mut seq = serializer.serialize_seq(Some(self.records.len()))?;
        for (record, &index) in self.records.iter().zip(&self.comparison.indices) {
            let s = &self.comparison.summaries[index];
            let legacy_sign = s.median.map(|v| {
                if v > 0.0 {
                    1
                } else if v < 0.0 {
                    -1
                } else {
                    0
                }
            });
            let mean = record.posterior_response_mean;
            let sign = if mean > 0.0 {
                1
            } else if mean < 0.0 {
                -1
            } else {
                0
            };
            let legacy = LegacySummary {
                lower: s.lower,
                mean: s.mean,
                median: s.median,
                mode: s.mode,
                posterior: (!self.comparison.certified).then_some(s.posterior.as_slice()),
                shared_posterior_index: self.comparison.certified.then_some(index),
                tail: s.tail,
                upper: s.upper,
            };
            seq.serialize_element(&ComparisonRecord {
                corrected_response: record,
                counts: record.counts,
                defined_differences: DEFINED_DIFFERENCES,
                descriptive_signs: Signs {
                    corrected_mean_sign: sign,
                    legacy_median_sign: legacy_sign,
                    signs_agree: legacy_sign.map(|s| s == sign),
                },
                direct_scientific_comparison: "unsupported",
                directional_interpretation_agrees: None,
                legacy_response: legacy,
                nt: &record.nt,
                reason: COMPARISON_REASON,
            })?;
        }
        seq.end()
    }
}
impl Prepared {
    fn data_view(&self) -> Clonotypes<'_, &[Clonotype]> {
        Clonotypes {
            header: &self.data,
            records: self.records.as_slice(),
        }
    }
    fn comparison_view(&self) -> Option<Clonotypes<'_, ComparisonRows<'_>>> {
        self.comparison.as_ref().map(|comparison| Clonotypes {
            header: &comparison.header,
            records: ComparisonRows {
                records: &self.records,
                comparison,
            },
        })
    }
    fn into_analysis(self) -> Result<Analysis> {
        let data =
            serde_json::to_value(self.data_view()).map_err(|e| err(State::NumericalFailure, e))?;
        let comparison = self
            .comparison_view()
            .map(serde_json::to_value)
            .transpose()
            .map_err(|e| err(State::NumericalFailure, e))?;
        Ok(Analysis {
            metadata: self.metadata,
            data,
            comparison,
        })
    }
}
fn float_text(value: f64) -> String {
    // Match Value::Number::to_string used by the accepted TSV writer,
    // including integral floats, scientific notation and signed zero.
    serde_json::Number::from_f64(value).map_or_else(|| "null".into(), |n| n.to_string())
}
const DEFINED_DIFFERENCES: &[&str] = &["legacy first-category versus corrected explicit-zero null conditioning", "legacy extra first-category factor versus corrected unseen mass in C2", "legacy one-sided extended f2 support versus bounded symmetric midpoint response prior", "legacy learned response hyperparameters versus caller-supplied corrected hyperparameters", "legacy median versus corrected posterior mean; no attribution of a clone difference to one cause"];
const COMPARISON_REASON: &str = "priors, support, fitted-versus-supplied response hyperparameters and summary definitions differ; descriptive signs are not biological calls";
fn response_error(e: ResponseError) -> AnalysisError {
    err(
        match e {
            ResponseError::PrecisionLoss | ResponseError::AllImpossible => State::UnresolvedResult,
            ResponseError::InvalidObservation => State::InvalidInput,
            _ => State::NumericalFailure,
        },
        e,
    )
}
fn validate(c: &Config) -> Result<()> {
    if c.operation != "fixed_parameter" || c.noise.model != "poisson" {
        return Err(err(State::UnsupportedOperation,"only fixed_parameter Poisson analysis is supported; no automatic fitting or model selection"));
    }
    if c.noise.parameter_order != ["exponent", "minimum_frequency"]
        || c.noise.parameters.iter().any(|v| !v.is_finite())
        || !(0.0 < c.noise.parameters[1] && c.noise.parameters[1] < 1.0)
        || !c.effect_delta.is_finite()
        || c.effect_delta < 0.0
        || c.noise.source_reference.trim().is_empty()
    {
        return Err(err(State::InvalidConfiguration,"require finite [exponent, linear minimum_frequency], 0<minimum<1, finite delta>=0, and a nonempty parameter source reference"));
    }
    if ![
        "user_supplied",
        "qualified_legacy_fit",
        "external",
        "experimental_corrected_fit",
    ]
    .contains(&c.noise.source.as_str())
    {
        return Err(err(State::InvalidConfiguration, "unknown parameter source"));
    }
    if c.noise.source != "qualified_legacy_fit" && c.noise.legacy_model_file.is_some() {
        return Err(err(
            State::InvalidConfiguration,
            "legacy_model_file requires qualified_legacy_fit source",
        ));
    }
    let p = &c.response_prior;
    if p.points < 2 || p.points > 1025 {
        return Err(err(
            State::InvalidConfiguration,
            "grid points must be 2..1025 under accepted prior contract",
        ));
    }
    Ok(())
}
fn source(c: &Config, base: &Path) -> Result<Value> {
    if c.noise.source != "qualified_legacy_fit" {
        return Ok(
            json!({"kind":c.noise.source,"reference":c.noise.source_reference,"hybrid":false,"estimation_qualification":"not assessed by corrected integration","statement":if c.noise.source=="experimental_corrected_fit"{"supplied values only; corrected estimator remains experimental and is not run here"}else{"caller-supplied values; estimation and biological validity are not qualified by this analysis"}}),
        );
    }
    let path = c.noise.legacy_model_file.as_ref().ok_or_else(|| {
        err(
            State::InvalidConfiguration,
            "qualified_legacy_fit requires a model.json source file",
        )
    })?;
    let path = base
        .join(path)
        .canonicalize()
        .map_err(|e| err(State::InvalidConfiguration, e))?;
    let bytes = fs::read(&path).map_err(|e| err(State::InvalidConfiguration, e))?;
    let model: Value =
        serde_json::from_slice(&bytes).map_err(|e| err(State::InvalidConfiguration, e))?;
    let x: [f64; 2] = serde_json::from_value(model["fit"]["evaluation"]["parameters"].clone())
        .map_err(|e| err(State::InvalidConfiguration, e))?;
    let residual = model["fit"]["evaluation"]["c2"]
        .as_f64()
        .unwrap_or(f64::NAN);
    if model["mode"] != "legacy-compatible"
        || model["upstream_revision"] != crate::legacy::REVISION
        || model["fit"]["success"] != true
        || x.iter().any(|v| !v.is_finite())
        || !residual.is_finite()
        || residual.abs() >= 1e-6
        || (x[0] - c.noise.parameters[0]).abs() > 1e-12
        || (x[1] - c.noise.parameters[1].log10()).abs() > 1e-12
    {
        return Err(err(
            State::InvalidConfiguration,
            "source must be a successful pinned legacy model matching declared exponent and cutoff",
        ));
    }
    Ok(
        json!({"kind":"qualified_legacy_fit","reference":c.noise.source_reference,"path":path,"sha256":crate::sha256(&bytes),"legacy_parameters_exponent_log10_minimum":x,"declared_parameters_exponent_linear_minimum":c.noise.parameters,"hybrid":true,"statement":"parameters estimated under legacy fitting contract; downstream inference uses corrected equations; hybrid analysis, not a fully corrected fitted workflow"}),
    )
}
/// In-memory integration used by file workflows and tests. Source files are
/// resolved relative to config_base. Rows must be unique nucleotide-key sorted.
pub fn analyze(
    rows: &[crate::MergedRow],
    c: &Config,
    config_base: &Path,
    compare: bool,
) -> Result<Analysis> {
    prepare_inner(rows, c, config_base, compare, false)?.into_analysis()
}
/// Explicit opt-in numerical certification, with the same frozen scientific configuration.
pub fn analyze_certified(
    rows: &[crate::MergedRow],
    c: &Config,
    config_base: &Path,
    compare: bool,
) -> Result<Analysis> {
    prepare_inner(rows, c, config_base, compare, true)?.into_analysis()
}
fn prepare_inner(
    rows: &[crate::MergedRow],
    c: &Config,
    config_base: &Path,
    compare: bool,
    certified: bool,
) -> Result<Prepared> {
    let mut profile = crate::profile::Profile::new("corrected-core");
    validate(c)?;
    let parameter_source = source(c, config_base)?;
    if rows.is_empty()
        || rows.windows(2).any(|w| w[0].nt >= w[1].nt)
        || rows
            .iter()
            .any(|r| r.nt.trim().is_empty() || (r.count1 == 0 && r.count2 == 0))
    {
        return Err(err(
            State::InvalidInput,
            "require sorted unique nonempty nucleotide keys and observed nonzero pairs",
        ));
    }
    let mut depth = [0_u64; 2];
    for row in rows {
        for (i, n) in [row.count1, row.count2].into_iter().enumerate() {
            depth[i] = depth[i]
                .checked_add(n)
                .ok_or_else(|| err(State::InvalidInput, "depth overflow"))?;
        }
    }
    if depth.contains(&0) || depth.iter().any(|d| *d > (1_u64 << 53)) {
        return Err(err(
            State::InvalidInput,
            "positive count-sum depths <=2^53 required",
        ));
    }
    let samples = depth.map(|d| SampleConfiguration {
        acquisition: Acquisition::Poisson,
        depth: d as f64,
    });
    let prior = build_symmetric_prior(
        Baseline {
            exponent: c.noise.parameters[0],
            minimum: c.noise.parameters[1],
        },
        c.second_baseline.map(|b| Baseline {
            exponent: b[0],
            minimum: b[1],
        }),
        ResponsePriorConfig {
            response_probability: c.response_prior.response_probability,
            scale: c.response_prior.scale,
            points: c.response_prior.points,
            stride: c.response_prior.stride,
            maximum_shift: c.response_prior.maximum_shift,
        },
    )
    .map_err(|e| {
        err(
            match e {
                crate::response_prior::PriorBuildError::InvalidConfiguration
                | crate::response_prior::PriorBuildError::IncompatibleBaselines => {
                    State::InvalidConfiguration
                }
                _ => State::NumericalFailure,
            },
            format!("symmetric prior: {e:?}"),
        )
    })?;
    let null = ConditionedModel::new(
        prior.grid.clone(),
        c.noise.parameters[0],
        samples[0],
        samples[1],
    )
    .map_err(|e| err(State::NumericalFailure, e))?;
    let table = PairTable::new(
        &rows
            .iter()
            .map(|r| ObservedCountPair {
                count1: r.count1,
                count2: r.count2,
                multiplicity: 1,
            })
            .collect::<Vec<_>>(),
    )
    .map_err(|e| err(State::InvalidInput, e))?;
    let obj = null
        .evaluate(&table)
        .map_err(|e| err(State::NumericalFailure, e))?;
    let posterior = ResponseModel::new(
        prior.prior.clone(),
        samples[0],
        samples[1],
        c.effect_delta,
        1.0,
    )
    .map_err(response_error)?;
    profile.mark("prior_null_objective_and_response_setup");
    // Scalar summaries suffice for file output. Do not retain each pair's
    // full posterior arrays; certification records remain independently kept.
    struct Summary {
        expanded: f64,
        contracted: f64,
        no_change: f64,
        mean_log_response: f64,
        credible_low: f64,
        credible_high: f64,
        log_evidence: f64,
        posterior_mass: f64,
        posterior_underflows: usize,
    }
    impl From<&crate::response::Posterior> for Summary {
        fn from(p: &crate::response::Posterior) -> Self {
            Self {
                expanded: p.expanded,
                contracted: p.contracted,
                no_change: p.no_change,
                mean_log_response: p.mean_log_response,
                credible_low: p.credible_low,
                credible_high: p.credible_high,
                log_evidence: p.log_evidence,
                posterior_mass: p.posterior_mass,
                posterior_underflows: p.posterior_underflows,
            }
        }
    }
    let mut records = Vec::new();
    let mut cache = BTreeMap::new();
    for row in rows {
        let key = (row.count1, row.count2);
        // Infer and diagnose in the original first-encounter row order. The
        // sorted cache changes neither arithmetic nor the first failing pair.
        if let std::collections::btree_map::Entry::Vacant(entry) = cache.entry(key) {
            let (p, certificate) = if certified {
                let (p, certificate) = posterior
                    .infer_certified_detailed(row.count1, row.count2)
                    .map_err(|failure| {
                    let mut error = response_error(failure.error);
                    if let Some(reason) = failure.reason {
                        error.message = format!("{}: {}", failure.error, reason.message());
                        error.certification_reason = Some(reason);
                    }
                    error
                })?;
                (p, Some(certificate))
            } else {
                (
                    posterior
                        .infer(row.count1, row.count2)
                        .map_err(response_error)?,
                    None,
                )
            };
            let pair = null
                .observed_pair(row.count1, row.count2)
                .map_err(|e| err(State::NumericalFailure, e))?;
            entry.insert((Summary::from(&p), pair, certificate));
        }
        let (p, pair, _) = &cache[&key];
        records.push(Clonotype {
            aa: row.aa.clone(),
            counts: [row.count1, row.count2],
            credible_interval_95: [p.credible_low, p.credible_high],
            log_evidence: p.log_evidence,
            negative_shift_probability: p.contracted,
            no_qualified_change_probability: p.no_change,
            nt: row.nt.clone(),
            null_pair: NullPair {
                log_conditioned: pair.conditioned.log_value(),
                log_posterior_mean_frequency: pair.log_posterior_mean_frequency,
                log_unconditioned: pair.unconditioned.log_value(),
            },
            positive_shift_probability: p.expanded,
            posterior_mass: p.posterior_mass,
            posterior_response_mean: p.mean_log_response,
            posterior_underflows: p.posterior_underflows,
            probability_within_delta: p.no_change,
            scientific_directional_call: None,
            state: "evaluated",
        });
    }
    profile.mark("posterior_null_pairs_and_record_assembly");
    let settings = json!({"noise":c.noise,"parameter_source":parameter_source,"response_prior":c.response_prior,"grid":{"points":c.response_prior.points,"minimum_frequency":c.noise.parameters[1],"maximum_frequency":1.0,"spacing":"uniform ln(f)","geometry":"u=ln(sqrt(f1*f2)), r=ln(f2/f1); pairs (i-k,i+k), diamond support"},"numerical":{"null_prior":"accepted analytically normalized power law, trapezoid masses discretely normalized","response_prior":"conditional midpoint slice normalization, exact mirrored masses, spike plus truncated discrete Laplace slab","effect_delta":c.effect_delta,"positive_event":"r>delta","negative_event":"r<-delta","no_qualified_change_probability_definition":"P(abs(r)<=delta), practical band; not probability that a decision rule failed","credible_interval":"accepted discrete equal-tail 95 percent","log_outputs_authoritative":true,"posterior_normalization_tolerance":3e-12,"prior_mass_tolerance":1e-12,"internal_call_threshold":1.0,"internal_call_exported":false,"precision_guard":"16*EPSILON*abs(acquisition log)<=5e-10"}});
    let mut metadata = json!({"analysis_mode":"corrected-fixed-parameter","noise_model":"poisson","settings":settings,"interpretation":"posterior probabilities conditional on supplied parameters, hyperparameters and discrete prior; null diagnostics refer to shared-frequency null model; no scientific directional calling rule","limitations":["no automatic corrected fitting or corrected response-hyperparameter learning","C1/C2 diagnostics only; no numeric C1 cutoff or fixed-parameter feasibility requirement","symmetric prior is bounded construction from accepted scout, not biological calibration","contraction probability is not a production contraction-selection rule","grid sensitivity/biological calibration not established for arbitrary inputs","supplied parameter uncertainty is not integrated","legacy-fit sources are hybrid, not fully corrected fitted analyses"],"version":crate::VERSION});
    if certified {
        metadata["numerical_certification"] = json!({"enabled":true,"requested":true,"omission_used":cache.values().any(|(_, _, c)| c.as_ref().is_some_and(|c| c.omitted_states > 0)),"failure_reason":null,"protocol_sha256":"6eb22de30367c3c97f584c71b416b64b6f1c47e6a6823951aa3e269dbcd3b196","mass_log2_budget":crate::certification::MASS_EXPONENT,"mean_log2_budget":crate::certification::MEAN_EXPONENT,"scope":"cumulative mathematical omission error; retained floating calculation separately qualified under unchanged accepted tolerances","certificate_index":"analysis.json:data.omission_certificates indexed by observed count pair","prior_support_changed":false});
    }
    let diagnostics = json!({"interpretation":"fixed-parameter shared-frequency null diagnostics, distinct from response joint-prior evidence","hidden_state":[0,0],"p00":obj.p00,"log_p00":obj.log_p00,"log_observation_probability":obj.log_observation_probability,"c1":obj.c1,"corrected_c2":obj.c2,"mean_nll":obj.mean_negative_log_likelihood,"total_nll":obj.negative_log_likelihood,"raw_prior_log_mass":obj.raw_prior_log_mass,"log_population_count":obj.log_population_count,"roundoff_corrections":obj.roundoff_corrections,"direct_probability_underflows":obj.direct_probability_underflows});
    let prior_data = json!({"raw_mass":prior.prior.raw_mass(),"state_count":prior.prior.states().len(),"states":prior.prior.states().iter().map(|s|json!({"frequency1":s.frequency1,"frequency2":s.frequency2,"mass":s.mass})).collect::<Vec<_>>(),"slices":prior.slices.iter().map(|s|json!({"shift":s.shift,"response":s.response,"one_signed_slice_mass":s.marginal_mass,"log_midpoint_normalizer":s.log_midpoint_normalizer,"midpoint_nodes":s.midpoint_nodes,"singleton":s.singleton})).collect::<Vec<_>>()});
    let mut comparison = None;
    if compare {
        let data = crate::legacy::Data::new(
            &rows
                .iter()
                .map(|r| [r.count1, r.count2])
                .collect::<Vec<_>>(),
        )
        .map_err(|e| err(State::UnsupportedOperation, e))?;
        let pars = [c.noise.parameters[0], c.noise.parameters[1].log10()];
        let legacy = crate::legacy::evaluate(&data, pars)
            .map_err(|e| err(State::LegacyReferenceFailure, e))?;
        let lr = crate::legacy::response(&data, pars, pars)
            .map_err(|e| err(State::LegacyReferenceFailure, e))?;
        let indices = rows
            .iter()
            .map(|row| {
                data.pairs
                    .iter()
                    .position(|p| p.counts == [row.count1, row.count2])
                    .unwrap()
            })
            .collect();
        let mut header = json!({"comparison_state":"evaluated","declared_parameters_exponent_linear_minimum":c.noise.parameters,"legacy_parameters_exponent_log10_minimum":pars,"whole_model":{"legacy":{"analysis_mode":"legacy-compatible","c1":legacy.c1,"c2":legacy.c2,"first_category":legacy.first_category,"conditioning_probability":legacy.first_probability,"explicit_zero_probability":legacy.explicit_zero_probability,"mean_nll":legacy.mean_nll,"grid_points":crate::legacy::POINTS,"response_grid":{"nodes":lr.shifts.len(),"minimum_shift":lr.shifts[0],"maximum_shift":lr.shifts[lr.shifts.len()-1]},"response_prior":{"definition":"legacy spike plus discrete Laplace; extended second-sample frequencies; nominal .1 reporting prior","learned_alpha":lr.alpha,"learned_scale":lr.scale},"limitations":"bounded compatibility; known conditioning/C2 and response-prior defects; no biological calibration"},"corrected":{"analysis_mode":"corrected-fixed-parameter","null_diagnostics":diagnostics,"response_prior_settings":c.response_prior,"response_state_count":prior.prior.states().len(),"minimum_response":prior.prior.states().iter().map(|s|s.frequency2.ln()-s.frequency1.ln()).fold(f64::INFINITY,f64::min),"maximum_response":prior.prior.states().iter().map(|s|s.frequency2.ln()-s.frequency1.ln()).fold(f64::NEG_INFINITY,f64::max),"limitations":metadata["limitations"]}},"clonotypes":null,"interpretation":"defined numerical/mathematical comparison only; differences do not establish biological truth"});
        if certified {
            let value = &mut header;
            value["shared_legacy_posteriors"]=json!(data.pairs.iter().zip(&lr.summaries).map(|(pair,s)|json!({"counts":pair.counts,"multiplicity":pair.multiplicity,"posterior":s.posterior})).collect::<Vec<_>>());
            value["certified_output_layout"]=json!("complete unchanged legacy posterior arrays stored once per count pair; per-clone summaries refer by shared_posterior_index; explicit certification mode only");
        }
        comparison = Some(PreparedComparison {
            header,
            summaries: lr.summaries,
            indices,
            certified,
        });
    }
    profile.mark("metadata_prior_and_optional_legacy_comparison");
    let mut result = Prepared {
        metadata,
        records,
        data: json!({"state":"evaluated","depths":depth,"observed_union":rows.len(),"null_diagnostics":diagnostics,"response_prior":null,"clonotypes":null}),
        comparison,
    };
    result.data["response_prior"] = prior_data;
    if certified {
        result.data["omission_certificates"] = json!(cache
            .iter()
            .map(|(&(a, b), (_, _, certificate))| { let certificate = certificate.as_ref().expect("certified cache entry"); json!({"counts":[a,b],"omission_used":certificate.omitted_states > 0,"interval_certification":if certificate.omitted_states > 0 {"passed"} else {"not_required_no_omission"},"certificate":certificate}) })
            .collect::<Vec<_>>());
    }
    profile.mark("result_value_assembly");
    Ok(result)
}
#[derive(Serialize)]
struct Envelope<'a, D: ?Sized> {
    data: &'a D,
    metadata: &'a Value,
}
// Hash below the buffering layer: serializer token writes are coalesced into
// bounded blocks before updating SHA256. An I/O failure cannot publish a stage.
struct StagedHash {
    file: fs::File,
    hasher: crate::output_hash::OutputHasher,
}
impl Write for StagedHash {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let written = self.file.write(bytes)?;
        self.hasher.update(&bytes[..written])?;
        Ok(written)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}
fn staged_error(context: &str, cause: impl fmt::Display) -> AnalysisError {
    err(State::NumericalFailure, format!("{context}: {cause}"))
}
fn staged_writer(stage: &Path, name: &str) -> Result<BufWriter<StagedHash>> {
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(stage.join(name))
        .map_err(|e| staged_error("create staged output", e))?;
    Ok(BufWriter::with_capacity(
        128 * 1024,
        StagedHash {
            file,
            hasher: crate::output_hash::OutputHasher::new(),
        },
    ))
}
fn finish_staged(mut writer: BufWriter<StagedHash>) -> Result<String> {
    writer
        .flush()
        .map_err(|e| staged_error("write staged output", e))?;
    let writer = writer
        .into_inner()
        .map_err(|e| staged_error("write staged output", e.error()))?;
    writer
        .file
        .sync_all()
        .map_err(|e| staged_error("write staged output", e))?;
    writer
        .hasher
        .finish()
        .map_err(|e| staged_error("hash staged output", e))
}
fn staged_json(stage: &Path, name: &str, value: &impl Serialize) -> Result<String> {
    let mut writer = staged_writer(stage, name)?;
    serde_json::to_writer_pretty(&mut writer, value).map_err(|e| {
        if e.is_io() {
            staged_error("write staged output", e)
        } else {
            err(State::NumericalFailure, e)
        }
    })?;
    writer
        .write_all(b"\n")
        .map_err(|e| staged_error("write staged output", e))?;
    finish_staged(writer)
}
/// Atomic publication only after every requested analysis completes. Every
/// artifact carries full mode/parameter/prior/settings/interpretation metadata.
pub fn run_files(
    first: &Path,
    second: &Path,
    config: &Path,
    output: &Path,
    compare: bool,
) -> Result<()> {
    run_files_inner(first, second, config, output, compare, false, None)
}
pub fn run_files_certified(
    first: &Path,
    second: &Path,
    config: &Path,
    output: &Path,
    compare: bool,
) -> Result<()> {
    run_files_inner(first, second, config, output, compare, true, None)
}
/// Internal integration retains streaming inference and its numerical certificates.
pub(crate) fn run_files_fitted(
    first: &Path,
    second: &Path,
    config: &Path,
    output: &Path,
    fitted: &Value,
) -> Result<()> {
    run_files_inner(first, second, config, output, false, true, Some(fitted))
}
fn run_files_inner(
    first: &Path,
    second: &Path,
    config: &Path,
    output: &Path,
    compare: bool,
    certified: bool,
    fitted: Option<&Value>,
) -> Result<()> {
    let mut profile = crate::profile::Profile::new("corrected-files");
    let a = crate::read_sample(first).map_err(|e| err(State::InvalidInput, e))?;
    profile.mark("sample1_read_hash_parse");
    let b = crate::read_sample(second).map_err(|e| err(State::InvalidInput, e))?;
    profile.mark("sample2_read_hash_parse");
    let rows = crate::merge(&a, &b).map_err(|e| err(State::InvalidInput, e))?;
    profile.mark("alignment");
    let config = config
        .canonicalize()
        .map_err(|e| err(State::InvalidConfiguration, e))?;
    let config_bytes = fs::read(&config).map_err(|e| err(State::InvalidConfiguration, e))?;
    let c: Config =
        serde_json::from_slice(&config_bytes).map_err(|e| err(State::InvalidConfiguration, e))?;
    let mut result = prepare_inner(&rows, &c, config.parent().unwrap(), compare, certified)?;
    result.metadata["input_provenance"] = json!([{"path":a.path,"sha256":a.sha256,"rows":a.rows,"depth":a.depth},{"path":b.path,"sha256":b.sha256,"rows":b.rows,"depth":b.depth}]);
    result.metadata["configuration_provenance"] = json!({"path":config,"sha256":crate::sha256(&config_bytes),"validation_policy":crate::POLICY});
    if let Some(fit) = fitted {
        result.metadata["analysis_mode"] = json!("corrected-auto-experimental");
        result.metadata["automatic_fit"] = fit.clone();
        result.metadata["settings"]["parameter_source"] = fit["parameter_source"].clone();
        result.metadata["configuration_provenance"]["path"] = fit["configuration_path"].clone();
        result.metadata["interpretation"] = json!("posterior probabilities conditional on experimental automatically fitted Poisson noise and supplied response hyperparameters; fitted parameter uncertainty is not integrated; no biological calibration or scientific calling rule");
        result.metadata["limitations"] = json!([
            "automatic noise fitting is EXPERIMENTAL; local constrained projection, no global optimum certificate",
            "training replicate roles are asserted by caller; biological changes can confound noise",
            "Poisson fit can be misspecified; split diagnostics do not qualify biological calibration",
            "response hyperparameters are supplied defaults or settings, not learned",
            "noise parameter uncertainty is not integrated into inference",
            "C1 remains a diagnostic without universal cutoff",
            "numerical omission certificates do not certify quadrature, estimation, or model validity",
            "no calibrated biological expansion/contraction calls"
        ]);
    }
    drop((a, b, rows));
    profile.mark("analysis_and_provenance");
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
        .canonicalize()
        .map_err(|e| err(State::InvalidConfiguration, e))?;
    let name = output.file_name().ok_or_else(|| {
        err(
            State::InvalidConfiguration,
            "output directory needs a new final name",
        )
    })?;
    let destination = parent.join(name);
    if fs::symlink_metadata(&destination).is_ok() {
        return Err(err(
            State::InvalidConfiguration,
            "refusing existing output directory",
        ));
    }
    let stage = tempfile::Builder::new()
        .prefix(".clonoshift-fixed-")
        .tempdir_in(parent)
        .map_err(|e| err(State::NumericalFailure, e))?;
    let mut hashes = BTreeMap::new();
    hashes.insert(
        "analysis.json",
        staged_json(
            stage.path(),
            "analysis.json",
            &Envelope {
                data: &result.data_view(),
                metadata: &result.metadata,
            },
        )?,
    );
    profile.mark("stream_analysis_json_and_sync");
    let mut table = staged_writer(stage.path(), "response.tsv")?;
    table
        .write_all(b"# metadata: ")
        .map_err(|e| staged_error("write staged output", e))?;
    serde_json::to_writer(&mut table, &result.metadata)
        .map_err(|e| staged_error("write staged output", e))?;
    table
        .write_all(b"\n")
        .map_err(|e| staged_error("write staged output", e))?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(table);
    writer
        .write_record([
            "ntCDR3",
            "count_1",
            "count_2",
            "positive_shift_probability",
            "negative_shift_probability",
            "probability_within_delta",
            "posterior_response_mean",
            "credible_low_95",
            "credible_high_95",
            "log_evidence",
            "state",
        ])
        .map_err(|e| err(State::NumericalFailure, e))?;
    for row in &result.records {
        let record = vec![
            row.nt.clone(),
            row.counts[0].to_string(),
            row.counts[1].to_string(),
            float_text(row.positive_shift_probability),
            float_text(row.negative_shift_probability),
            float_text(row.probability_within_delta),
            float_text(row.posterior_response_mean),
            float_text(row.credible_interval_95[0]),
            float_text(row.credible_interval_95[1]),
            float_text(row.log_evidence),
            "evaluated".into(),
        ];
        writer
            .write_record(record)
            .map_err(|e| err(State::NumericalFailure, e))?;
    }
    let table = writer
        .into_inner()
        .map_err(|e| err(State::NumericalFailure, e))?;
    hashes.insert("response.tsv", finish_staged(table)?);
    profile.mark("serialize_response_tsv");
    let mut workflow_metadata = result.metadata.clone();
    if let Some(comparison) = result.comparison {
        let mut meta = result.metadata.clone();
        meta["analysis_mode"] = json!("legacy-versus-corrected");
        meta["mode_definitions"] = json!({"legacy":"legacy equations, response hyperparameters learned under legacy path","corrected":"supplied fixed noise parameters and response hyperparameters, accepted corrected kernels"});
        meta["legacy_settings"] = comparison.header["whole_model"]["legacy"].clone();
        workflow_metadata = meta.clone();
        if certified {
            hashes.insert("comparison-summary.json",staged_json(stage.path(),"comparison-summary.json",&json!({"metadata":meta,"data":{
                "whole_model":comparison.header["whole_model"],
                "observed_union":Some(result.records.len()),
                "shared_legacy_posterior_groups":comparison.header["shared_legacy_posteriors"].as_array().map(Vec::len),
                "complete_rows_and_posteriors":"comparison.json",
                "interpretation":comparison.header["interpretation"],
                "layout":comparison.header["certified_output_layout"]
            }}))?);
        }
        hashes.insert(
            "comparison.json",
            staged_json(
                stage.path(),
                "comparison.json",
                &Envelope {
                    data: &Clonotypes {
                        header: &comparison.header,
                        records: ComparisonRows {
                            records: &result.records,
                            comparison: &comparison,
                        },
                    },
                    metadata: &meta,
                },
            )?,
        );
    }
    profile.mark("serialize_optional_comparison");
    staged_json(
        stage.path(),
        "provenance.json",
        &json!({"metadata":workflow_metadata,"data":{"state":"evaluated","workflow":if fitted.is_some(){"corrected-auto-experimental"}else if compare{"compare-fixed"}else{"corrected-fixed"},"outputs_sha256":hashes,"atomic_publication":"no-replace directory rename; no partial successful output","noise_parameter_estimation_run":fitted.is_some(),"legacy_response_hyperparameter_learning_run":compare,"corrected_response_hyperparameter_learning_run":false,"scientific_calls_emitted":false}}),
    )?;
    profile.mark("stream_provenance_json_and_sync");
    fs::File::open(stage.path())
        .and_then(|f| f.sync_all())
        .map_err(|e| err(State::NumericalFailure, e))?;
    crate::publish(stage.path(), &destination).map_err(|e| err(State::NumericalFailure, e))?;
    profile.mark("synchronized_atomic_publication");
    Ok(())
}
pub fn error_json(e: &AnalysisError) -> String {
    let mut value = json!({"state":e.state,"message":e.message});
    if let Some(reason) = e.certification_reason {
        value["classification"] = json!("PrecisionLoss");
        value["reason_code"] = json!(reason);
    }
    value.to_string()
}

#[cfg(test)]
mod serialization_tests {
    use super::*;
    #[test]
    fn typed_file_views_preserve_original_keys_float_bytes_and_legacy_layouts() {
        let row = Clonotype {
            aa: "Q\"Ω".into(),
            counts: [0, 7],
            credible_interval_95: [-0.0, 1e300],
            log_evidence: -1e300,
            negative_shift_probability: 1e-300,
            no_qualified_change_probability: 0.0,
            nt: "A\"Ω".into(),
            null_pair: NullPair {
                log_conditioned: -0.0,
                log_posterior_mean_frequency: -1.5,
                log_unconditioned: -2.0,
            },
            positive_shift_probability: 1.0,
            posterior_mass: 1.0,
            posterior_response_mean: 1.5,
            posterior_underflows: 3,
            probability_within_delta: 0.0,
            scientific_directional_call: None,
            state: "evaluated",
        };
        let original = json!({"AACDR3":"Q\"Ω","ntCDR3":"A\"Ω","counts":[0,7],
            "credible_interval_95":[-0.0,1e300],"log_evidence":-1e300,
            "negative_shift_probability":1e-300,"no_qualified_change_probability":0.0,
            "null_pair":{"log_conditioned":-0.0,"log_posterior_mean_frequency":-1.5,"log_unconditioned":-2.0},
            "positive_shift_probability":1.0,"posterior_mass":1.0,"posterior_response_mean":1.5,
            "posterior_underflows":3,"probability_within_delta":0.0,
            "scientific_directional_call":null,"state":"evaluated"});
        assert_eq!(
            serde_json::to_vec_pretty(&row).unwrap(),
            serde_json::to_vec_pretty(&original).unwrap()
        );
        let posterior = [-0.0, 1e-300, 1.0];
        for certified in [false, true] {
            let summary = LegacySummary {
                lower: None,
                mean: 1.5,
                median: Some(-0.0),
                mode: 0.0,
                posterior: (!certified).then_some(posterior.as_slice()),
                shared_posterior_index: certified.then_some(7),
                tail: 1e-300,
                upper: Some(1e300),
            };
            let mut expected = json!({"lower":null,"mean":1.5,"median":-0.0,"mode":0.0,"tail":1e-300,"upper":1e300});
            if certified {
                expected["shared_posterior_index"] = json!(7);
            } else {
                expected["posterior"] = json!(posterior);
            }
            assert_eq!(
                serde_json::to_vec_pretty(&summary).unwrap(),
                serde_json::to_vec_pretty(&expected).unwrap()
            );
        }
        assert_eq!(float_text(-0.0), "-0.0");
        assert_eq!(float_text(1.0), "1.0");
    }
}
