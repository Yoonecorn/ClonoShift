# ClonoShift User Manual

ClonoShift 1.0.0.

Qualified platform: macOS on Apple Silicon.

License: GPL-3.0-or-later.

This guide covers the standard source-build workflow.

## 1. What ClonoShift does

ClonoShift compares repertoire counts from two samples. A repertoire is a collection of immune-receptor sequences and their counts. You define which receptor unit each count represents.

Extended and Core analyze the same qualified input. Here, qualified input means counts that meet the software's requirements. It does not mean that the biological interpretation has been established.

The two models use different assumptions and summaries. A model is a mathematical description of how the observations can arise. ClonoShift reports the two results separately. It does not select a winner, average the models, or make a final biological call.

Follow the [installation procedure](#2-install-clonoshift), then [run the example](#3-run-the-example). Use the example to check the workflow before you prepare study data.

## 2. Install ClonoShift

Use Terminal on your Mac. You need the downloaded or extracted ClonoShift 1.0.0 repository, Rust, and Cargo. Cargo is the tool that builds the executable. The verified local toolchain uses Rust and Cargo 1.98.0.

You also need the existing local dependency cache. A dependency cache stores the software packages required for the build. Cargo must already be configured to use that cache in your installation.

The repository root contains `prepare_config.rb`. This supplied Ruby helper is required to create the analysis configuration. Ruby and its JSON library must be available. The helper was checked with Ruby 2.6.10.

### Build and verify

Open Terminal in the downloaded or extracted repository root. Copy the complete block below to build the executable and check its version.

```bash
export CARGO_TARGET_DIR="$PWD/target"
cargo build --release --locked --offline
target/release/clonoshift --version
```

The expected version is:

```text
clonoshift 1.0.0
```

The executable is `target/release/clonoshift`. Keep Terminal in the repository root for the example. Paths such as `fixtures/synthetic/counts.tsv` start from this folder.

The build uses the recorded dependency versions and the local cache. If Cargo reports a missing offline package, ask your installation maintainer to supply or configure the required cache.

The `CARGO_TARGET_DIR` line puts the executable in this source folder. Keep that line even if your computer has another Cargo build location configured.

A backslash at the end of a command line means that the command continues on the next line. Copy the full command block. Do not put a space after a backslash.

## 3. Run the example

The bundled example uses synthetic counts and arbitrary identifiers. It checks the workflow. It does not establish biological findings.

Start in the repository root. The example creates `fit`, `hybrid.json`, and `results`. These names must not already exist. Keep earlier outputs and use a fresh working copy if you repeat the exercise.

### Fit the example noise model

A noise model describes variation between measurements. Run:

```bash
target/release/clonoshift legacy-poisson \
  --sample1 fixtures/synthetic/replicate1.tsv \
  --sample2 fixtures/synthetic/replicate2.tsv --output-dir fit
```

Expect `legacy-compatible output: fit`. The folder contains the fitted `model.json` used by the next steps.

### Create the configuration

A configuration file stores the model settings. Run the required helper:

```bash
ruby prepare_config.rb fit/model.json hybrid.json
```

Expect a message that begins `Created hybrid.json`. The helper uses the successful fit and copies the supplied example prior settings. A prior describes the responses allowed by a model before it uses the response counts.

The example prior settings are illustrative. They need scientific review before use with real data. The helper does not choose a suitable biological prior for your study.

### Run Extended and Core

Run both models with the bundled merged count file:

```bash
target/release/clonoshift analyze \
  --counts fixtures/synthetic/counts.tsv --output results \
  --config hybrid.json --model1 fit/model.json --model2 fit/model.json
```

Both models run by default. The example uses the same synthetic fit for both model inputs. A real study needs noise models suitable for its samples and assay.

A successful `analyze` command prints no success message. Check the saved statuses instead of assuming that silence means both models completed.

### Check the statuses and find the results

Open the run record with:

```bash
cat results/comparison/run-metadata.json
```

JSON is a structured text format. Find `extended`, then its `status`. Find `core`, then its `status`. Both must say `COMPLETED` for this example. If a model says `REFUSED`, read its `reason` and use [Common problems](#7-common-problems).

The four main result files are:

| File inside `results` | What to read |
|---|---|
| `comparison/run-metadata.json` | Model statuses, reasons, and run settings |
| `extended/response.tsv` | Extended results for each analyzed unit |
| `core/response.tsv` | Core results for each analyzed unit |
| `comparison/summary.tsv` | Counts of descriptive sign agreement, disagreement, and unavailability |

Read the example summary:

```bash
cat results/comparison/summary.tsv
```

Expect agreement `5`, disagreement `1`, and unavailable `0`. These are comparisons of numerical signs. They are not counts of confirmed biological findings.

## 4. Prepare your data

ClonoShift matches the identifiers you supply. It does not establish their biological meaning, reconstruct chain pairing, or infer clonal families.

**Caution:** Use the same biological-unit definition in both samples. Do not mix nucleotide clonotypes, amino-acid clonotypes, exact heavy/light pairs, or clonal families.

### Sample files

Prepare one tab-separated values file, or TSV file, for each sample. A tab separates one field from the next. Save plain text rather than a spreadsheet workbook.

A minimal sample file is:

```text
ntCDR3	count
NT001	12
NT002	4
```

Use one row for each defined unit. `ntCDR3` is the identifier used to match that unit between samples. Its name refers to nucleotide complementarity-determining region 3, a variable receptor region. The parser matches the supplied identifier; it does not establish its biological meaning.

Counts must be whole numbers of zero or more. Use the original integer counts, not normalized frequencies. Identifiers must be unique within each file. Each sample must have a total count greater than zero.

The optional `AACDR3` column holds an amino-acid annotation. If both sample files give an annotation for the same identifier, the nonempty annotations must agree. You can omit this column when you do not have an annotation.

Keep count meaning consistent. Reads, molecule counts, and cell counts are different measurements. Do not combine them as if they were interchangeable.

### Merged analysis file

The standard validation command creates the merged file from the two sample files. A merged file looks like this:

```text
ntCDR3	AACDR3	count_1	count_2
NT001	CASSF	12	18
NT002	CASRF	0	4
```

Use these four columns in this order. Separate fields with tabs. The `AACDR3` column must be present in this merged format, but its values can be empty.

`count_1` belongs to sample 1. `count_2` belongs to sample 2. Keep that sample order throughout the study, including the noise-model inputs.

Do not use duplicate identifiers. Do not analyze a row with zero counts in both samples. A unit absent from one sample can have zero in that sample and a positive count in the other.

Use the validation step below to align identifiers between files. Do not align rows by their positions in a spreadsheet. Check the reported sample totals against your original data before analysis.

## 5. Run your analysis

In the standard workflow, Core fits measurement-noise parameters from suitable replicate data. Parameters are numerical values used by a model. Extended uses parameters from a suitable successful Core fit and holds them fixed during analysis.

Technical replicates are repeated measurements used to fit measurement noise. Sample 1 and sample 2 define the biological comparison you want to analyze. Do not treat a biological comparison as technical replication without scientific justification.

Start in the built repository root from Section 2. Create a separate study folder so the example outputs remain intact:

```bash
mkdir study
```

Put `sample1.tsv`, `sample2.tsv`, `replicate1.tsv`, and `replicate2.tsv` in this `study` folder before you continue. The first two files contain the response samples. The other two files contain suitable technical replicates for the noise fit.

Run all commands below from the repository root. Paths under `study/` refer to the study folder. The output names `study/prepared`, `study/fit`, `study/study.json`, and `study/results` must be new.

### Validate and merge the response samples

```bash
target/release/clonoshift validate --sample1 study/sample1.tsv \
  --sample2 study/sample2.tsv --output-dir study/prepared
```

Expect a message beginning `validated`. Check the two reported sample totals. The merged analysis input is `study/prepared/merged.tsv`.

If validation fails, correct the source data before proceeding. Do not change counts merely to make the program accept them.

### Fit the noise model

```bash
target/release/clonoshift legacy-poisson \
  --sample1 study/replicate1.tsv --sample2 study/replicate2.tsv --output-dir study/fit
```

Expect `legacy-compatible output: study/fit`. Continue only with a successful fit. The fitted model is `study/fit/model.json`.

### Create the configuration

Use the helper with the successful model and the supplied example template:

```bash
ruby prepare_config.rb study/fit/model.json study/study.json
```

Expect a message beginning `Created study/study.json`. The helper copies the fit values into a compatible configuration. Extended holds these fitted parameters fixed during response analysis.

The copied prior settings are illustrative and require scientific review before real-data use. This includes `effect_delta`, the size of the band treated as practically close to zero. Obtain suitable settings for your assay and study before running an analysis for interpretation.

The helper creates the file. It does not establish that the replicate design, noise model, or response prior is appropriate.

### Run both models

```bash
target/release/clonoshift analyze --counts study/prepared/merged.tsv \
  --output study/results --config study/study.json \
  --model1 study/fit/model.json --model2 study/fit/model.json
```

**Important:** Do not reuse one fitted noise model for both samples unless that calibration is scientifically suitable for both.

This command uses one fit twice. If the samples need different calibrations, supply their suitable model files in sample order. Extended's fitted parameters must match a suitable successful Core fit. That fitted model must match one of the supplied Core models.

### Check both statuses

```bash
cat study/results/comparison/run-metadata.json
```

Read `extended.status` and `core.status`. Use a model's estimates only when its status is `COMPLETED`. Read the reason for every other status before interpreting the run.

Keep the complete `results` folder. It includes provenance: records of the inputs, settings, and file identities used for the run. Keep the source sample files, fit, and configuration with your study records.

## 6. Read the results

Open `comparison/run-metadata.json` inside the output folder first. This record tells you which model results are available.

| Status | Meaning |
|---|---|
| `COMPLETED` | The model produced a usable result under its stated assumptions. |
| `REFUSED` | The model withheld a result because a required condition was not met. |
| `FAILED` | Processing or calculation failed. |
| `UNAVAILABLE` | The model was not requested or could not be compared. |

A command can finish without reporting a command failure while a model is `REFUSED`. Always read both statuses. A completed result still depends on suitable data, calibration, and model settings.

The response describes the change in a unit's frequency from sample 1 to sample 2. Frequency means its share of the repertoire. Positive response means higher inferred frequency in sample 2. Negative response means lower inferred frequency. Zero means equal frequencies.

### Extended

Read `extended/response.tsv`. Its first line contains metadata, or information about the run. The next line gives the column names. Each following row represents one analyzed unit.

A posterior describes model probabilities after the model uses the counts. Read these main fields:

| Field | Plain meaning |
|---|---|
| `posterior_response_mean` | The probability-weighted average response over Extended's response grid, its finite set of evaluated states. |
| `credible_low_95` | The lower end of the model's 95% credible interval. |
| `credible_high_95` | The upper end of that interval. |
| `probability_within_delta` | Probability that the response lies inside the supplied practical-equivalence band around zero. |
| `positive_shift_probability` | Probability that the response is above the positive edge of that band. |
| `negative_shift_probability` | Probability that the response is below the negative edge of that band. |

A credible interval contains a stated share of posterior probability, subject to the model's discrete grid. An interval that includes zero is compatible with zero response under that model.

Practical equivalence means change inside a chosen band around zero. Its probability is not a final statement that nothing happened biologically. All three region probabilities depend on the supplied band and model settings.

### Core

Read `core/response.tsv`. Its columns describe Core's own response distribution. Read these fields:

| Field | Plain meaning |
|---|---|
| `median` | Core's discrete posterior median, a middle location of its response distribution. |
| `lower` | The lower reported 95% interval bound. |
| `upper` | The upper reported 95% interval bound. |
| `tail_P_s_le_zero` | Probability of a response at or below zero, including probability at exactly zero. |

The Core tail is not a pure contraction probability. A reported value of `unavailable` for an interval or median is not zero.

### Read the models together

Both models use the same qualified counts. Their numerical values can differ because they use different assumptions and summaries. Extended reports a mean; Core reports a discrete median.

Agreement can support robustness across the two formulations. A zero median and a small nonzero mean are not automatically contradictory. They can describe distributions with substantial probability near zero.

Do not average the results or select whichever result better matches an expected biological story. Read each model under its own definitions. A valid result from one model can remain usable if the other model refuses.

The comparison summary describes agreement between the signs of Extended's mean and Core's median. It does not supply a third biological result. A small nonzero mean and a zero median can be listed as disagreement.

## 7. Common problems

### Command not found

Terminal cannot find the executable at the supplied location. Return to the stated starting folder and use the exact executable path shown for that step.

### Input file not found

`No such file or directory` can mean the input path is wrong. Check that the named file is in the current study folder, with the same spelling and letter case.

### Required column missing

The file does not have the required header. Restore the sample or merged header shown in [Prepare your data](#4-prepare-your-data), including the merged column order.

### Invalid count

`invalid count` means a count does not meet the integer requirements. Replace it with the correct original whole-number count; do not substitute a normalized frequency.

### Duplicate identifier

`duplicate nucleotide key` means an identifier occurs more than once in a sample file. Resolve the repeated records using your original data and a documented rule for the biological unit.

### Output folder already exists

`refusing existing output` means that the destination is already present, even if empty. Choose a new output name and keep the earlier result.

### One model returned REFUSED

A required model condition was not met. Read that model's `reason` in `comparison/run-metadata.json` and check the named data, model, or configuration requirement.

### Extended and Core differ

The models can give different summaries of the same counts. Check sample order and settings, then read both results using their separate definitions in [Read the results](#6-read-the-results).

## 8. Important limits and quick reference

ClonoShift does not provide a consensus biological call. Posterior probabilities are not p-values. Model agreement does not prove biological correctness.

Count data alone cannot identify whether extra variation is technical or biological. Calibration must match the assay and study design. Calibration means checking a model against suitable measurements or a stated reference law.

A refusal is not a zero result. Do not turn an unresolved estimate into a finding of no change.

B-cell receptor, or BCR, binding evidence applies to the assayed exact heavy/light pair. Do not assign that evidence automatically to related unassayed sequences. A related clonal-family member is not the same assayed receptor pair.

Use this checklist with the standard commands in [Run your analysis](#5-run-your-analysis):

1. Verify the version.
2. Validate the sample files.
3. Fit the noise model.
4. Create the configuration.
5. Run both models.
6. Check model status.
7. Read Extended and Core separately.
8. Keep the complete output folder and provenance.
