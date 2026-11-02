// SPDX-License-Identifier: GPL-3.0-or-later
//! Experimental replicate training followed by certified Poisson inference.
//! Training roles are an explicit caller assertion, distinct from target samples.
use crate::analysis::{AnalysisError, Config, Noise, Prior, State};
use crate::observation::{ObservedCountPair, PairTable};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufReader, Read, Write};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Settings {
    pub response_prior: Prior,
    pub effect_delta: f64,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            response_prior: Prior {
                response_probability: 0.2,
                scale: 1.0,
                points: 129,
                stride: 1,
                maximum_shift: 24,
            },
            effect_delta: 0.3,
        }
    }
}
fn error(state: State, message: impl ToString) -> AnalysisError {
    AnalysisError {
        state,
        message: message.to_string(),
        certification_reason: None,
    }
}
fn write_json(path: &Path, value: &impl Serialize) -> Result<(), AnalysisError> {
    let mut file = fs::File::create(path).map_err(|e| error(State::NumericalFailure, e))?;
    serde_json::to_writer(&mut file, value).map_err(|e| error(State::NumericalFailure, e))?;
    file.write_all(b"\n")
        .and_then(|()| file.sync_all())
        .map_err(|e| error(State::NumericalFailure, e))
}
fn file_hash(path: &Path) -> Result<String, AnalysisError> {
    use sha2::{Digest, Sha256};
    let mut stream =
        BufReader::new(fs::File::open(path).map_err(|e| error(State::NumericalFailure, e))?);
    let mut hash = Sha256::new();
    let mut buffer = [0_u8; 65536];
    loop {
        let n = stream
            .read(&mut buffer)
            .map_err(|e| error(State::NumericalFailure, e))?;
        if n == 0 {
            break;
        }
        hash.update(&buffer[..n]);
    }
    Ok(format!("{:x}", hash.finalize()))
}
/// Runs only after explicit CLI experimental opt-in. A refused fit publishes a
/// diagnostic bundle with no inference and returns UnresolvedResult. Success
/// means an experimental calculation completed, never production qualification.
pub fn run_files(
    replicates: [&Path; 2],
    samples: [&Path; 2],
    settings: &Settings,
    output: &Path,
) -> Result<(), AnalysisError> {
    // Validate output before expensive input/fit work; never overwrite evidence.
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
        .canonicalize()
        .map_err(|e| error(State::InvalidConfiguration, e))?;
    let name = output.file_name().ok_or_else(|| {
        error(
            State::InvalidConfiguration,
            "new final output name required",
        )
    })?;
    let destination = parent.join(name);
    if fs::symlink_metadata(&destination).is_ok() {
        return Err(error(
            State::InvalidConfiguration,
            "refusing existing output directory",
        ));
    }
    let a = crate::read_sample(replicates[0]).map_err(|e| error(State::InvalidInput, e))?;
    let b = crate::read_sample(replicates[1]).map_err(|e| error(State::InvalidInput, e))?;
    if a.sha256 == b.sha256 {
        return Err(error(
            State::InvalidInput,
            "identical replicate bytes do not supply independent noise evidence",
        ));
    }
    let depths = [a.depth, b.depth];
    if depths.contains(&0) || depths.iter().any(|n| *n > (1_u64 << 53)) {
        return Err(error(
            State::InvalidInput,
            "positive replicate depths <=2^53 required",
        ));
    }
    let rows = crate::merge(&a, &b).map_err(|e| error(State::InvalidInput, e))?;
    let provenance = json!([
        {"path":a.path,"sha256":a.sha256,"rows":a.rows,"depth":a.depth},
        {"path":b.path,"sha256":b.sha256,"rows":b.rows,"depth":b.depth}
    ]);
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
    .map_err(|e| error(State::InvalidInput, e))?;
    drop((a, b, rows));
    let report = crate::autofit::fit(
        &table,
        depths.map(|n| n as f64),
        &crate::autofit::AutoFitConfig::default(),
    );
    drop(table);
    let report_value =
        serde_json::to_value(&report).map_err(|e| error(State::NumericalFailure, e))?;
    let stage = tempfile::Builder::new()
        .prefix(".clonoshift-auto-")
        .tempdir_in(&parent)
        .map_err(|e| error(State::NumericalFailure, e))?;
    let fit_path = destination.join("fit-report.json");
    let report_envelope = json!({"classification":"EXPERIMENTAL","training_inputs":provenance,
        "replicate_role":"caller-asserted same repertoire replicates; independence and absence of biological change not verified from count data",
        "fit":report_value,"requested_target_inputs":samples,"response_settings":settings,"response_hyperparameters_learned":false,
        "biological_calibration":false,"parameter_uncertainty_integrated":false});
    write_json(&stage.path().join("fit-report.json"), &report_envelope)?;
    let config_path = stage.path().join("inferred-fixed-config.json");
    let candidate = &report_value["candidate"];
    let state = if report_value["outcome"] == "Candidate" {
        if let (Some(exponent), Some(minimum)) = (
            candidate["exponent"].as_f64(),
            candidate["minimum"].as_f64(),
        ) {
            let config = Config {
                noise: Noise {
                    model: "poisson".into(),
                    parameter_order: ["exponent".into(), "minimum_frequency".into()],
                    parameters: [exponent, minimum],
                    source: "experimental_corrected_fit".into(),
                    source_reference: fit_path.to_string_lossy().into(),
                    legacy_model_file: None,
                },
                response_prior: settings.response_prior.clone(),
                effect_delta: settings.effect_delta,
                second_baseline: None,
                operation: "fixed_parameter".into(),
            };
            write_json(&config_path, &config)?;
            let context = json!({"classification":"EXPERIMENTAL","configuration_path":destination.join("inferred-fixed-config.json"),
                "fit_report":fit_path,"fit_report_sha256":file_hash(&stage.path().join("fit-report.json"))?,
                "parameter_source":{"kind":"experimental_corrected_fit","reference":fit_path,"hybrid":false,
                    "estimation_qualification":"Experimental C2 constrained Poisson projection; see complete fit report",
                    "noise_parameter_estimation_run":true,"training_inputs":provenance},
                "parameter_uncertainty_integrated":false,"response_hyperparameters_learned":false});
            crate::analysis::run_files_fitted(
                samples[0],
                samples[1],
                &config_path,
                &stage.path().join("inference"),
                &context,
            )
        } else {
            Err(error(
                State::UnresolvedResult,
                "fitter reported Candidate without finite parameters",
            ))
        }
    } else {
        Err(error(
            State::UnresolvedResult,
            format!("automatic fit refused: {}", report_value["reason"]),
        ))
    };
    if let Err(ref e) = state {
        write_json(
            &stage.path().join("failure.json"),
            &json!({"state":e.state,"message":e.message,
            "certification_reason":e.certification_reason.map(|r|format!("{r:?}")),"inference_published":false}),
        )?;
    }
    let mut hashes = BTreeMap::<String, String>::new();
    for name in [
        "fit-report.json",
        "inferred-fixed-config.json",
        "failure.json",
    ] {
        let path = stage.path().join(name);
        if path.is_file() {
            hashes.insert(name.into(), file_hash(&path)?);
        }
    }
    let inference = stage.path().join("inference");
    if inference.exists() {
        for entry in fs::read_dir(inference).map_err(|e| error(State::NumericalFailure, e))? {
            let entry = entry.map_err(|e| error(State::NumericalFailure, e))?;
            if entry.path().is_file() {
                hashes.insert(
                    format!("inference/{}", entry.file_name().to_string_lossy()),
                    file_hash(&entry.path())?,
                );
            }
        }
    }
    write_json(
        &stage.path().join("provenance.json"),
        &json!({"workflow":"corrected-auto-experimental","classification":"EXPERIMENTAL",
        "state":if state.is_ok(){"evaluated"}else{"refused"},"noise_parameter_estimation_run":true,
        "response_hyperparameter_learning_run":false,"scientific_calls_emitted":false,"outputs_sha256":hashes,
        "atomic_publication":"no-replace directory rename; failed calculations retain diagnostic bundle only"}),
    )?;
    fs::File::open(stage.path())
        .and_then(|f| f.sync_all())
        .map_err(|e| error(State::NumericalFailure, e))?;
    crate::publish(stage.path(), &destination).map_err(|e| error(State::NumericalFailure, e))?;
    state
}
