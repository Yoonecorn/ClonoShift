// SPDX-License-Identifier: GPL-3.0-or-later
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};

fn execute(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_clonoshift"))
        .args(args)
        .output()
        .unwrap()
}
fn path(p: &Path) -> &str {
    p.to_str().unwrap()
}
#[test]
fn dual_workflow_replays_and_preserves_other_result_on_refusal() {
    let work = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fixture = root.join("fixtures/synthetic");
    let fit = work.path().join("fit");
    let z = execute(&[
        "legacy-poisson",
        "--sample1",
        path(&fixture.join("replicate1.tsv")),
        "--sample2",
        path(&fixture.join("replicate2.tsv")),
        "--output-dir",
        path(&fit),
    ]);
    assert!(z.status.success(), "{}", String::from_utf8_lossy(&z.stderr));
    let model = fit.join("model.json");
    let v: Value = serde_json::from_slice(&fs::read(&model).unwrap()).unwrap();
    let params = v["fit"]["evaluation"]["parameters"].as_array().unwrap();
    let config = work.path().join("config.json");
    let c = json!({"noise":{"model":"poisson","parameter_order":["exponent","minimum_frequency"],"parameters":[params[0].as_f64().unwrap(),10.0_f64.powf(params[1].as_f64().unwrap())],"source":"qualified_legacy_fit","source_reference":"synthetic fixture","legacy_model_file":path(&model)},"response_prior":{"response_probability":0.2,"scale":1.0,"points":33,"stride":1,"maximum_shift":6},"effect_delta":0.3});
    fs::write(&config, serde_json::to_vec(&c).unwrap()).unwrap();
    let counts = fixture.join("counts.tsv");
    let run = |out: &Path, cfg: &Path| {
        execute(&[
            "analyze",
            "--counts",
            path(&counts),
            "--output",
            path(out),
            "--config",
            path(cfg),
            "--model1",
            path(&model),
            "--model2",
            path(&model),
        ])
    };
    let a = work.path().join("a");
    let b = work.path().join("b");
    for p in [&a, &b] {
        let z = run(p, &config);
        assert!(z.status.success(), "{}", String::from_utf8_lossy(&z.stderr));
    }
    for p in [
        "extended/response.tsv",
        "core/response.tsv",
        "comparison/agreement.tsv",
        "comparison/disagreement.tsv",
        "comparison/run-metadata.json",
    ] {
        assert_eq!(
            fs::read(a.join(p)).unwrap(),
            fs::read(b.join(p)).unwrap(),
            "{p}"
        );
    }
    for model_name in ["extended", "core"] {
        let actual = fs::read_to_string(a.join(model_name).join("response.tsv")).unwrap();
        let actual = actual
            .lines()
            .filter(|line| !line.starts_with("# metadata:"))
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        let expected =
            fs::read_to_string(fixture.join("expected").join(format!("{model_name}.tsv"))).unwrap();
        assert_eq!(actual, expected, "{model_name} numerical table");
    }
    assert!(fs::read_to_string(a.join("comparison/run-metadata.json"))
        .unwrap()
        .starts_with("{\"extended\":"));
    assert!(!run(&a, &config).status.success());
    let bad = work.path().join("bad.json");
    fs::write(&bad, br#"{"noise":{"source":"user_supplied"}}"#).unwrap();
    let refused = work.path().join("refused");
    let z = run(&refused, &bad);
    assert!(z.status.success());
    let m: Value =
        serde_json::from_slice(&fs::read(refused.join("comparison/run-metadata.json")).unwrap())
            .unwrap();
    assert_eq!(m["extended"]["status"], "REFUSED");
    assert_eq!(m["core"]["status"], "COMPLETED");
    assert!(refused.join("core/response.tsv").exists());
    let missing = work.path().join("missing-model.json");
    let failed = work.path().join("failed");
    let z = execute(&[
        "analyze",
        "--counts",
        path(&counts),
        "--output",
        path(&failed),
        "--config",
        path(&config),
        "--model1",
        path(&model),
        "--model2",
        path(&missing),
    ]);
    assert!(!z.status.success());
    let m: Value =
        serde_json::from_slice(&fs::read(failed.join("comparison/run-metadata.json")).unwrap())
            .unwrap();
    assert_eq!(m["extended"]["status"], "COMPLETED");
    assert_eq!(m["core"]["status"], "FAILED");
    assert!(failed.join("extended/response.tsv").exists());
    let bad_model = work.path().join("invalid-model.json");
    fs::write(&bad_model, b"{}").unwrap();
    let refused_core = work.path().join("refused-core");
    let z = execute(&[
        "analyze",
        "--counts",
        path(&counts),
        "--output",
        path(&refused_core),
        "--config",
        path(&config),
        "--model1",
        path(&model),
        "--model2",
        path(&bad_model),
    ]);
    assert!(z.status.success());
    let m: Value = serde_json::from_slice(
        &fs::read(refused_core.join("comparison/run-metadata.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(m["extended"]["status"], "COMPLETED");
    assert_eq!(m["core"]["status"], "REFUSED");
    let hash_key = work.path().join("hash-key.tsv");
    fs::write(
        &hash_key,
        fs::read_to_string(&counts)
            .unwrap()
            .replace("NT000", "#NT000"),
    )
    .unwrap();
    let hash_out = work.path().join("hash-key-out");
    let z = execute(&[
        "analyze",
        "--counts",
        path(&hash_key),
        "--output",
        path(&hash_out),
        "--config",
        path(&config),
        "--model1",
        path(&model),
        "--model2",
        path(&model),
    ]);
    assert!(z.status.success());
    let compare = fs::read_to_string(hash_out.join("comparison/agreement.tsv")).unwrap()
        + &fs::read_to_string(hash_out.join("comparison/disagreement.tsv")).unwrap()
        + &fs::read_to_string(hash_out.join("comparison/unavailable.tsv")).unwrap();
    assert!(compare.contains("#NT000\t"));
    let invalid = work.path().join("invalid.tsv");
    fs::write(&invalid, b"ntCDR3\tAACDR3\tcount_1\tcount_2\nx\tA\t-1\t2\n").unwrap();
    let output = work.path().join("invalid-output");
    let z = execute(&[
        "analyze",
        "--counts",
        path(&invalid),
        "--output",
        path(&output),
        "--config",
        path(&config),
        "--model1",
        path(&model),
        "--model2",
        path(&model),
    ]);
    assert!(!z.status.success());
    assert!(!output.exists());
}
