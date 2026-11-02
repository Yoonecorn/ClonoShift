// SPDX-License-Identifier: GPL-3.0-or-later
mod analyze;
mod automatic;
mod fixed;
use std::path::PathBuf;

const LEGACY_USAGE: &str = "Legacy compatibility: clonoshift legacy-poisson --sample1 FILE --sample2 FILE --output-dir NEW_DIRECTORY\nResponse: clonoshift legacy-response --sample1 FILE --sample2 FILE --model1 MODEL_JSON --model2 MODEL_JSON --output-dir NEW_DIRECTORY\nKnown pinned NoisET defects are isolated in these commands. See docs/METHODS.md.\n";
const USAGE: &str = "ClonoShift — validated repertoire ingestion\n\nUsage: clonoshift validate --sample1 FILE --sample2 FILE --output-dir NEW_DIRECTORY\n\nTSV headers: ntCDR3,count[,AACDR3] or N. Seq. CDR3,Clone count[,AA. Seq. CDR3].\nOutput parent must exist; destination must not exist. No files are overwritten.\n";
fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args_os()
        .skip(1)
        .map(|value| {
            value
                .into_string()
                .map_err(|value| format!("command arguments must be UTF-8: {value:?}"))
        })
        .collect::<Result<_, _>>()?;
    if args.as_slice() == ["--help"] || args.as_slice() == ["-h"] {
        print!(
            "{}{USAGE}{LEGACY_USAGE}{}{}",
            analyze::USAGE,
            fixed::USAGE,
            automatic::USAGE
        );
        return Ok(());
    }
    if args.as_slice() == ["--version"] {
        println!("clonoshift {}", clonoshift::VERSION);
        return Ok(());
    }
    if args.first().is_some_and(|s| s == "analyze") {
        return analyze::run(&args);
    }
    if args.first().is_some_and(|s| s == "corrected-auto") {
        return automatic::run(&args);
    }
    if args
        .first()
        .is_some_and(|s| s == "corrected-fixed" || s == "compare-fixed")
    {
        return fixed::run(&args);
    }
    if args
        .first()
        .is_some_and(|s| s == "legacy-poisson" || s == "legacy-response")
    {
        return legacy(&args);
    }
    if args.first().map(String::as_str) != Some("validate") {
        return Err(USAGE.to_owned());
    }
    let mut paths: [Option<PathBuf>; 3] = [None, None, None];
    let mut index = 1;
    while index < args.len() {
        if args[index] == "--help" {
            print!("{USAGE}{LEGACY_USAGE}");
            return Ok(());
        }
        let slot = match args[index].as_str() {
            "--sample1" => 0,
            "--sample2" => 1,
            "--output-dir" => 2,
            other => return Err(format!("unknown option {other:?}\n{USAGE}")),
        };
        if paths[slot].is_some() {
            return Err(format!("duplicate option {}", args[index]));
        }
        let value = args
            .get(index + 1)
            .filter(|v| !v.starts_with("--"))
            .ok_or_else(|| format!("missing value for {}", args[index]))?;
        paths[slot] = Some(value.into());
        index += 2;
    }
    let [Some(first), Some(second), Some(output)] = paths else {
        return Err(format!("all three options are required\n{USAGE}"));
    };
    let result = clonoshift::validate(&first, &second, &output).map_err(|e| e.to_string())?;
    println!(
        "validated {} union clonotypes; depths {} / {}; output {}",
        result.union_clonotype_count,
        result.inputs[0].depth,
        result.inputs[1].depth,
        result.output_path.display()
    );
    Ok(())
}
fn legacy(args: &[String]) -> Result<(), String> {
    let response = args[0] == "legacy-response";
    let mut paths: [Option<PathBuf>; 5] = [None, None, None, None, None];
    if args.len() == 2 && args[1] == "--help" {
        print!("{LEGACY_USAGE}");
        return Ok(());
    }
    let mut i = 1;
    while i < args.len() {
        let slot = match args[i].as_str() {
            "--sample1" => 0,
            "--sample2" => 1,
            "--output-dir" => 2,
            "--model1" if response => 3,
            "--model2" if response => 4,
            _ => {
                return Err(format!(
                    "unsupported legacy option {}\n{LEGACY_USAGE}",
                    args[i]
                ))
            }
        };
        if paths[slot].is_some() {
            return Err("duplicate legacy option".into());
        }
        let value = args
            .get(i + 1)
            .filter(|v| !v.starts_with("--"))
            .ok_or("missing legacy option value")?;
        paths[slot] = Some(value.into());
        i += 2;
    }
    if paths[..if response { 5 } else { 3 }]
        .iter()
        .any(Option::is_none)
    {
        return Err(format!("missing required legacy options\n{LEGACY_USAGE}"));
    }
    clonoshift::legacy::run_files(
        paths[0].as_ref().unwrap(),
        paths[1].as_ref().unwrap(),
        paths[2].as_ref().unwrap(),
        if response {
            Some([
                paths[3].as_ref().unwrap().as_path(),
                paths[4].as_ref().unwrap().as_path(),
            ])
        } else {
            None
        },
    )?;
    println!(
        "legacy-compatible output: {}",
        paths[2].as_ref().unwrap().display()
    );
    Ok(())
}
fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("clonoshift: {message}");
            std::process::ExitCode::from(2)
        }
    }
}
