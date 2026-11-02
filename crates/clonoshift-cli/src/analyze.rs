// SPDX-License-Identifier: GPL-3.0-or-later
//! Dual-model file orchestration. Scientific calculations remain in the existing modules.
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs,
    io::Write,
    path::{Path, PathBuf},
};

pub const USAGE: &str = "Analyze: clonoshift analyze --counts MERGED_TSV --output NEW_DIRECTORY [--model both|extended|core] --config HYBRID_JSON --model1 MODEL_JSON --model2 MODEL_JSON\nOptional --certify-negligible-states enables the qualified omission certificate. The counts file has ntCDR3, AACDR3, count_1, count_2 columns. Config is required for Extended; both models are required for Core.\n";

fn bytes(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}
fn save(path: &Path, data: &[u8]) -> Result<(), String> {
    let mut f = fs::File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    f.write_all(data)
        .and_then(|_| f.sync_all())
        .map_err(|e| format!("{}: {e}", path.display()))
}
fn table(path: &Path) -> Result<BTreeMap<String, Vec<String>>, String> {
    let s = String::from_utf8(bytes(path)?).map_err(|e| e.to_string())?;
    let mut lines = s.lines().filter(|l| !l.starts_with("# metadata:"));
    let header = lines.next().ok_or("missing table header")?;
    let cols = header.split('\t').count();
    let mut out = BTreeMap::new();
    for line in lines {
        let row: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if row.len() != cols || out.insert(row[0].clone(), row).is_some() {
            return Err("malformed or duplicate output row".into());
        }
    }
    Ok(out)
}
fn compare(
    root: &Path,
    ext: &str,
    core: &str,
    ex_reason: &str,
    co_reason: &str,
) -> Result<(), String> {
    let mut agree = String::from("ntCDR3\textended_mean\tcore_median\n");
    let mut disagree = agree.clone();
    let mut unavailable = String::from("ntCDR3\treason\n");
    let mut a = 0;
    let mut d = 0;
    let mut u = 0;
    if ext == "COMPLETED" && core == "COMPLETED" {
        let left = table(&root.join("extended/response.tsv"))?;
        let right = table(&root.join("core/response.tsv"))?;
        if left.keys().ne(right.keys()) {
            return Err("model result keys differ".into());
        }
        for (key, x) in &left {
            let y = &right[key];
            let em = x[6].parse::<f64>().map_err(|e| e.to_string())?;
            match y[7].parse::<f64>() {
                Ok(cm) if em.is_finite() && cm.is_finite() => {
                    let row = format!("{key}\t{}\t{}\n", x[6], y[7]);
                    if (em > 0.0) == (cm > 0.0) && (em < 0.0) == (cm < 0.0) {
                        agree.push_str(&row);
                        a += 1;
                    } else {
                        disagree.push_str(&row);
                        d += 1;
                    }
                }
                _ => {
                    unavailable.push_str(&format!("{key}\tcore median unavailable\n"));
                    u += 1;
                }
            }
        }
    } else {
        unavailable.push_str(&format!(
            "*\textended={ext}: {}; core={core}: {}\n",
            serde_json::to_string(ex_reason).unwrap(),
            serde_json::to_string(co_reason).unwrap()
        ));
        u += 1;
    }
    let p = root.join("comparison");
    fs::create_dir(&p).map_err(|e| e.to_string())?;
    save(&p.join("agreement.tsv"), agree.as_bytes())?;
    save(&p.join("disagreement.tsv"), disagree.as_bytes())?;
    save(&p.join("unavailable.tsv"), unavailable.as_bytes())?;
    save(
        &p.join("summary.tsv"),
        format!("metric\tvalue\nagreement\t{a}\ndisagreement\t{d}\nunavailable\t{u}\n").as_bytes(),
    )?;
    Ok(())
}

fn normalize_provenance(root: &Path, scratch: &Path, model: &str) -> Result<(), String> {
    let dir = root.join(model);
    for name in if model == "extended" {
        vec!["analysis.json", "response.tsv", "provenance.json"]
    } else {
        vec!["provenance.json"]
    } {
        let path = dir.join(name);
        let mut text = String::from_utf8(bytes(&path)?).map_err(|e| e.to_string())?;
        for n in ["sample1.tsv", "sample2.tsv"] {
            let original = scratch.join(n).canonicalize().map_err(|e| e.to_string())?;
            text = text.replace(
                original.to_str().ok_or("non-UTF-8 temporary path")?,
                &format!("derived/{n}"),
            );
        }
        save(&path, text.as_bytes())?;
    }
    if model == "extended" {
        let path = dir.join("provenance.json");
        let mut value: Value = serde_json::from_slice(&bytes(&path)?).map_err(|e| e.to_string())?;
        for n in ["analysis.json", "response.tsv"] {
            value["data"]["outputs_sha256"][n] = json!(clonoshift::sha256(&bytes(&dir.join(n))?));
        }
        let mut encoded = serde_json::to_vec_pretty(&value).map_err(|e| e.to_string())?;
        encoded.push(b'\n');
        save(&path, &encoded)?;
    }
    Ok(())
}

fn core_run(
    first: &Path,
    second: &Path,
    output: &Path,
    models: &[PathBuf; 2],
) -> Result<(), (&'static str, String)> {
    for model in models {
        let raw = bytes(model).map_err(|e| ("FAILED", e))?;
        let value: Value = serde_json::from_slice(&raw).map_err(|e| ("REFUSED", e.to_string()))?;
        if value["mode"] != "legacy-compatible"
            || value["upstream_revision"] != clonoshift::legacy::REVISION
            || value["fit"]["success"] != true
        {
            return Err((
                "REFUSED",
                "response requires successful legacy-compatible model.json files".into(),
            ));
        }
        let _: [f64; 2] = serde_json::from_value(value["fit"]["evaluation"]["parameters"].clone())
            .map_err(|e| ("REFUSED", e.to_string()))?;
    }
    clonoshift::legacy::run_files(first, second, output, Some([&models[0], &models[1]])).map_err(
        |e| {
            let refused = [
                "unrepresentable legacy prior support",
                "legacy arithmetic outside supported finite range",
                "invalid legacy response prior",
                "unrepresentable legacy response prior",
                "invalid legacy posterior arrays",
                "zero legacy posterior evidence",
                "unsupported legacy response shift stride",
                "nonfinite legacy response surface",
            ]
            .contains(&e.as_str());
            (if refused { "REFUSED" } else { "FAILED" }, e)
        },
    )
}

pub fn run(args: &[String]) -> Result<(), String> {
    if args.len() == 2 && args[1] == "--help" {
        print!("{USAGE}");
        return Ok(());
    }
    let mut opts = BTreeMap::new();
    let mut i = 1;
    let mut certified = false;
    while i < args.len() {
        if args[i] == "--certify-negligible-states" {
            if certified {
                return Err("duplicate --certify-negligible-states".into());
            }
            certified = true;
            i += 1;
            continue;
        }
        let k = args[i].as_str();
        if ![
            "--counts", "--output", "--model", "--config", "--model1", "--model2",
        ]
        .contains(&k)
        {
            return Err(format!("unknown option {k}\n{USAGE}"));
        }
        let v = args
            .get(i + 1)
            .filter(|v| !v.starts_with("--"))
            .ok_or_else(|| format!("missing value for {k}"))?;
        if opts.insert(k, v.as_str()).is_some() {
            return Err(format!("duplicate option {k}"));
        }
        i += 2;
    }
    let get = |k| {
        opts.get(k)
            .copied()
            .ok_or_else(|| format!("missing {k}\n{USAGE}"))
    };
    let input = PathBuf::from(get("--counts")?);
    let output = PathBuf::from(get("--output")?);
    let mode = opts.get("--model").copied().unwrap_or("both");
    if !["both", "extended", "core"].contains(&mode) {
        return Err("invalid --model".into());
    }
    let config = if mode != "core" {
        Some(PathBuf::from(get("--config")?))
    } else {
        None
    };
    let models = if mode != "extended" {
        Some([
            PathBuf::from(get("--model1")?),
            PathBuf::from(get("--model2")?),
        ])
    } else {
        None
    };
    let input_bytes = bytes(&input)?;
    let text = std::str::from_utf8(&input_bytes).map_err(|e| e.to_string())?;
    let mut lines = text.lines();
    if lines.next() != Some("ntCDR3\tAACDR3\tcount_1\tcount_2") {
        return Err("counts header must be ntCDR3, AACDR3, count_1, count_2 in that order".into());
    }
    let mut first = String::from("ntCDR3\tAACDR3\tcount\n");
    let mut second = first.clone();
    for line in lines {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 4 || fields.iter().any(|x| x.contains('"') || x.contains('\r')) {
            return Err("invalid unquoted counts record".into());
        }
        first.push_str(&format!("{}\t{}\t{}\n", fields[0], fields[1], fields[2]));
        second.push_str(&format!("{}\t{}\t{}\n", fields[0], fields[1], fields[3]));
    }
    let one =
        clonoshift::parse_sample(first.as_bytes(), "sample1".into()).map_err(|e| e.to_string())?;
    let two =
        clonoshift::parse_sample(second.as_bytes(), "sample2".into()).map_err(|e| e.to_string())?;
    let rows = clonoshift::merge(&one, &two).map_err(|e| e.to_string())?;
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
        .canonicalize()
        .map_err(|e| e.to_string())?;
    let name = output.file_name().ok_or("invalid output name")?;
    let destination = parent.join(name);
    if fs::symlink_metadata(&destination).is_ok() {
        return Err("refusing existing output".into());
    }
    let stage = tempfile::Builder::new()
        .prefix(".clonoshift-analyze-")
        .tempdir_in(&parent)
        .map_err(|e| e.to_string())?;
    // Keep staging on any subsequent error so completed model artifacts survive.
    let root_path = stage.keep();
    let root = root_path.as_path();
    let outcome = (|| -> Result<(), String> {
        let scratch = tempfile::tempdir().map_err(|e| e.to_string())?;
        let p1 = scratch.path().join("sample1.tsv");
        let p2 = scratch.path().join("sample2.tsv");
        save(&p1, first.as_bytes())?;
        save(&p2, second.as_bytes())?;
        let mut ex = "UNAVAILABLE";
        let mut co = "UNAVAILABLE";
        let mut ex_reason = "not requested".to_string();
        let mut co_reason = "not requested".to_string();
        let mut config_hash = None;
        let mut model_hashes: Vec<String> = Vec::new();
        if let Some(m) = &models {
            for path in m {
                if let Ok(b) = bytes(path) {
                    model_hashes.push(clonoshift::sha256(&b));
                }
            }
        }
        if let Some(c) = &config {
            let attempt = std::panic::catch_unwind(|| -> Result<(), (String, String)> {
                let refuse = |e: String| ("REFUSED".to_string(), e);
                let cfg_bytes = bytes(c).map_err(refuse)?;
                let cfg: Value =
                    serde_json::from_slice(&cfg_bytes).map_err(|e| refuse(e.to_string()))?;
                if cfg["noise"]["source"] != "qualified_legacy_fit" {
                    return Err(refuse("Extended requires qualified_legacy_fit".into()));
                }
                let model_path = cfg["noise"]["legacy_model_file"]
                    .as_str()
                    .ok_or_else(|| refuse("missing legacy_model_file".into()))?;
                let cfg_parent = c
                    .canonicalize()
                    .map_err(|e| refuse(e.to_string()))?
                    .parent()
                    .ok_or_else(|| refuse("config has no parent".into()))?
                    .to_path_buf();
                let digest =
                    clonoshift::sha256(&bytes(&cfg_parent.join(model_path)).map_err(refuse)?);
                if let Some(m) = &models {
                    let hashes = m
                        .iter()
                        .filter_map(|p| bytes(p).ok().map(|b| clonoshift::sha256(&b)))
                        .collect::<Vec<_>>();
                    if !hashes.contains(&digest) {
                        return Err(refuse("Extended fit must match one Core model".into()));
                    }
                }
                let result = if certified {
                    clonoshift::analysis::run_files_certified(
                        &p1,
                        &p2,
                        c,
                        &root.join("extended"),
                        false,
                    )
                } else {
                    clonoshift::analysis::run_files(&p1, &p2, c, &root.join("extended"), false)
                };
                result.map_err(|e| {
                    let status = if matches!(e.state, clonoshift::analysis::State::NumericalFailure)
                    {
                        "FAILED"
                    } else {
                        "REFUSED"
                    };
                    (status.to_string(), clonoshift::analysis::error_json(&e))
                })
            });
            config_hash = bytes(c).ok().map(|b| clonoshift::sha256(&b));
            match attempt {
                Ok(Ok(())) => {
                    ex = "COMPLETED";
                    ex_reason.clear();
                }
                Ok(Err((status, reason))) => {
                    ex = if status == "FAILED" {
                        "FAILED"
                    } else {
                        "REFUSED"
                    };
                    ex_reason = reason;
                }
                Err(_) => {
                    ex = "FAILED";
                    ex_reason = "unexpected panic during Extended".into();
                }
            }
            if ex != "COMPLETED" {
                fs::create_dir_all(root.join("extended")).map_err(|e| e.to_string())?;
            }
        } else {
            fs::create_dir(root.join("extended")).map_err(|e| e.to_string())?;
        }
        if let Some(m) = &models {
            let pairs: Vec<_> = rows.iter().map(|r| [r.count1, r.count2]).collect();
            match clonoshift::legacy::Data::new(&pairs) {
                Err(e) => {
                    co = "REFUSED";
                    co_reason = e;
                    fs::create_dir(root.join("core")).map_err(|e| e.to_string())?;
                }
                Ok(_) => {
                    match std::panic::catch_unwind(|| core_run(&p1, &p2, &root.join("core"), m)) {
                        Ok(Ok(())) => {
                            co = "COMPLETED";
                            co_reason.clear();
                        }
                        Ok(Err((state, e))) => {
                            co = state;
                            co_reason = e;
                            fs::create_dir_all(root.join("core")).map_err(|e| e.to_string())?;
                        }
                        Err(_) => {
                            co = "FAILED";
                            co_reason = "unexpected panic during Core".into();
                            fs::create_dir_all(root.join("core")).map_err(|e| e.to_string())?;
                        }
                    }
                }
            }
        } else {
            fs::create_dir(root.join("core")).map_err(|e| e.to_string())?;
        }
        if ex == "COMPLETED" {
            if let Err(e) = normalize_provenance(root, scratch.path(), "extended") {
                ex = "FAILED";
                ex_reason = format!("provenance normalization: {e}");
            }
        }
        if co == "COMPLETED" {
            if let Err(e) = normalize_provenance(root, scratch.path(), "core") {
                co = "FAILED";
                co_reason = format!("provenance normalization: {e}");
            }
        }
        let comparison_error = compare(root, ex, co, &ex_reason, &co_reason).err();
        if let Some(reason) = &comparison_error {
            fs::create_dir_all(root.join("comparison")).map_err(|e| e.to_string())?;
            for name in ["summary.tsv", "agreement.tsv", "disagreement.tsv"] {
                let path = root.join("comparison").join(name);
                if !path.exists() {
                    save(&path, b"status\nFAILED\n")?;
                }
            }
            save(
                &root.join("comparison/unavailable.tsv"),
                format!(
                    "ntCDR3\treason\n*\t{}\n",
                    serde_json::to_string(reason).unwrap()
                )
                .as_bytes(),
            )?;
        }
        let config_parameters = config
            .as_ref()
            .and_then(|p| bytes(p).ok())
            .and_then(|b| serde_json::from_slice::<Value>(&b).ok());
        let core_parameters: Vec<Value> = models
            .as_ref()
            .map(|m| {
                m.iter()
                    .map(|p| {
                        bytes(p)
                            .ok()
                            .and_then(|b| serde_json::from_slice::<Value>(&b).ok())
                            .map(|v| v["fit"]["evaluation"]["parameters"].clone())
                            .unwrap_or(Value::Null)
                    })
                    .collect()
            })
            .unwrap_or_default();
        let meta = json!({"version":clonoshift::VERSION,"input_sha256":clonoshift::sha256(&input_bytes),"derived_sample_sha256":[clonoshift::sha256(first.as_bytes()),clonoshift::sha256(second.as_bytes())],"model":mode,"certify_negligible_states":certified,"comparison_error":comparison_error,"extended":{"status":ex,"reason":ex_reason,"config_sha256":config_hash,"config":config_parameters},"core":{"status":co,"reason":co_reason,"model_sha256":model_hashes,"parameters":core_parameters},"comparison":"descriptive signs of Extended mean and Core median, different estimands; agreement is not biological truth","warnings":["No biological winner, calibrated directional calls or biological FDR is inferred."]});
        let mut remaining = meta
            .as_object()
            .ok_or("internal metadata is not an object")?
            .clone();
        let extended = remaining
            .remove("extended")
            .ok_or("missing Extended metadata")?;
        let core = remaining.remove("core").ok_or("missing Core metadata")?;
        let rest = serde_json::to_string(&remaining).map_err(|e| e.to_string())?;
        let metadata = format!(
            "{{\"extended\":{},\"core\":{},{}\n",
            serde_json::to_string(&extended).map_err(|e| e.to_string())?,
            serde_json::to_string(&core).map_err(|e| e.to_string())?,
            &rest[1..]
        );
        save(
            &root.join("comparison/run-metadata.json"),
            metadata.as_bytes(),
        )?;
        for dir in ["extended", "core", "comparison"] {
            fs::File::open(root.join(dir))
                .and_then(|f| f.sync_all())
                .map_err(|e| e.to_string())?;
        }
        fs::File::open(root)
            .and_then(|f| f.sync_all())
            .map_err(|e| e.to_string())?;
        clonoshift::publish(root, &destination).map_err(|e| e.to_string())?;
        if ex == "FAILED" || co == "FAILED" || comparison_error.is_some() {
            return Err("analysis failed; completed artifacts were published".into());
        }
        Ok(())
    })();
    outcome.map_err(|e| {
        if root_path.exists() {
            format!("{e}; preserved staging at {}", root_path.display())
        } else {
            e
        }
    })
}
