// SPDX-License-Identifier: GPL-3.0-or-later
use clonoshift::analysis::{AnalysisError, State};
use clonoshift::automatic::Settings;
use std::path::PathBuf;
pub const USAGE: &str = "Experimental: clonoshift corrected-auto --experimental --replicate1 FILE --replicate2 FILE --sample1 FILE --sample2 FILE --output-dir NEW_DIRECTORY [--settings JSON]\nReplicates train Poisson noise; sample1/2 receive certified corrected inference. Settings supply response_prior and effect_delta only. Fitting is EXPERIMENTAL; parameter uncertainty and biological calibration unavailable. Refusals retain diagnostic bundles only.\n";
pub fn run(args: &[String]) -> Result<(), String> {
    let action = || -> Result<(), AnalysisError> {
        let invalid = |message: &str| AnalysisError {
            state: State::InvalidConfiguration,
            message: message.into(),
            certification_reason: None,
        };
        if args.len() == 2 && args[1] == "--help" {
            print!("{USAGE}");
            return Ok(());
        }
        let mut paths: [Option<PathBuf>; 6] = [None, None, None, None, None, None];
        let mut experimental = false;
        let mut i = 1;
        while i < args.len() {
            if args[i] == "--experimental" {
                if experimental {
                    return Err(invalid("duplicate experimental option"));
                }
                experimental = true;
                i += 1;
                continue;
            }
            let slot = match args[i].as_str() {
                "--replicate1" => 0,
                "--replicate2" => 1,
                "--sample1" => 2,
                "--sample2" => 3,
                "--output-dir" => 4,
                "--settings" => 5,
                _ => return Err(invalid("unsupported automatic workflow option")),
            };
            if paths[slot].is_some() {
                return Err(invalid("duplicate automatic workflow option"));
            }
            let value = args
                .get(i + 1)
                .filter(|s| !s.starts_with("--"))
                .ok_or_else(|| invalid("missing option value"))?;
            paths[slot] = Some(value.into());
            i += 2;
        }
        if !experimental {
            return Err(invalid(
                "Corrected automatic fitting is EXPERIMENTAL; explicit --experimental required",
            ));
        }
        let [Some(a), Some(b), Some(x), Some(y), Some(out), settings_path] = paths else {
            return Err(invalid(
                "replicate1, replicate2, sample1, sample2 and output-dir required",
            ));
        };
        let settings = if let Some(path) = settings_path {
            let bytes = std::fs::read(path).map_err(|e| invalid(&e.to_string()))?;
            serde_json::from_slice::<Settings>(&bytes).map_err(|e| invalid(&e.to_string()))?
        } else {
            Settings::default()
        };
        clonoshift::automatic::run_files([&a, &b], [&x, &y], &settings, &out)?;
        println!("experimental automatic corrected output: {}", out.display());
        Ok(())
    };
    action().map_err(|e| clonoshift::analysis::error_json(&e))
}
