// SPDX-License-Identifier: GPL-3.0-or-later
use clonoshift::analysis::{AnalysisError, State};
use std::path::PathBuf;
pub const USAGE:&str="Corrected fixed parameters: clonoshift corrected-fixed --sample1 FILE --sample2 FILE --config JSON --output-dir NEW_DIRECTORY\nComparison: clonoshift compare-fixed (same options). No automatic corrected fitting or biological directional calls. See docs/METHODS.md.\nOptional --certify-negligible-states enables numerical omission certification; see docs/METHODS.md.\n";
pub fn run(args: &[String]) -> Result<(), String> {
    let action = || -> Result<(), AnalysisError> {
        if args.len() == 2 && args[1] == "--help" {
            print!("{USAGE}");
            return Ok(());
        }
        let mut paths: [Option<PathBuf>; 4] = [None, None, None, None];
        let mut i = 1;
        let mut certified = false;
        let invalid = |m: &str| AnalysisError {
            state: State::InvalidConfiguration,
            message: m.into(),
            certification_reason: None,
        };
        while i < args.len() {
            if args[i] == "--certify-negligible-states" {
                if certified {
                    return Err(invalid("duplicate certification option"));
                }
                certified = true;
                i += 1;
                continue;
            }
            let slot = match args[i].as_str() {
                "--sample1" => 0,
                "--sample2" => 1,
                "--config" => 2,
                "--output-dir" => 3,
                _ => return Err(invalid("unsupported corrected workflow option")),
            };
            if paths[slot].is_some() {
                return Err(invalid("duplicate corrected workflow option"));
            }
            let value = args
                .get(i + 1)
                .filter(|s| !s.starts_with("--"))
                .ok_or_else(|| invalid("missing option value"))?;
            paths[slot] = Some(value.into());
            i += 2;
        }
        let [Some(a), Some(b), Some(config), Some(out)] = paths else {
            return Err(invalid(
                "sample1, sample2, config and output-dir are required",
            ));
        };
        if certified {
            clonoshift::analysis::run_files_certified(
                &a,
                &b,
                &config,
                &out,
                args[0] == "compare-fixed",
            )?;
        } else {
            clonoshift::analysis::run_files(&a, &b, &config, &out, args[0] == "compare-fixed")?;
        }
        println!("corrected fixed-parameter output: {}", out.display());
        Ok(())
    };
    action().map_err(serde_json_error)
}
fn serde_json_error(e: AnalysisError) -> String {
    // The CLI crate already depends on the core; serialize classified failures there.
    clonoshift::analysis::error_json(&e)
}
