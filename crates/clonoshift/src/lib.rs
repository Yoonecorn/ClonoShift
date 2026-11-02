// SPDX-License-Identifier: GPL-3.0-or-later
//! Validated repertoire ingestion, independently qualified acquisition kernels,
//! observation conditioning, fixed noise objectives and explicit-prior posteriors.
//! Bounded fitting and midpoint-prior scouts are separate from accepted kernels.
pub mod acquisition;
pub mod analysis;
pub mod autofit;
pub mod automatic;
pub mod certification;
pub mod fitting;
pub mod legacy;
pub mod numerics;
pub mod observation;
mod output_hash;
mod profile;
pub mod response;
pub mod response_objective;
pub mod response_prior;

use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const POLICY: &str = "exact UTF-8 nucleotide key; no normalization; empty/whitespace-only keys rejected; ASCII decimal u64 counts; checked u64 depths; duplicates rejected; conflicting nonempty AA rejected; absent counts zero; exclusivity by key presence; frequencies ignored; UTF-8 lexicographic output order";

#[derive(Debug)]
pub struct Error(String);
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}
impl std::error::Error for Error {}
pub type Result<T> = std::result::Result<T, Error>;
fn error(message: impl Into<String>) -> Error {
    Error(message.into())
}
fn context(message: &str, cause: impl fmt::Display) -> Error {
    error(format!("{message}: {cause}"))
}
pub fn sha256(bytes: &[u8]) -> String {
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    if let Some(digest) = apple_sha256(bytes) {
        const HEX: &[u8; 16] = b"0123456789abcdef";
        let mut encoded = String::with_capacity(64);
        for byte in digest {
            encoded.push(HEX[(byte >> 4) as usize] as char);
            encoded.push(HEX[(byte & 15) as usize] as char);
        }
        return encoded;
    }
    format!("{:x}", Sha256::digest(bytes))
}

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
fn apple_sha256(bytes: &[u8]) -> Option<[u8; 32]> {
    // Public CommonDigest.h: CC_LONG is uint32_t. Larger slices keep the
    // portable implementation; never truncate a length at the FFI boundary.
    let length = u32::try_from(bytes.len()).ok()?;
    #[link(name = "System")]
    unsafe extern "C" {
        fn CC_SHA256(data: *const libc::c_void, len: u32, digest: *mut u8) -> *mut u8;
    }
    let mut digest = [0_u8; 32];
    // SAFETY: The immutable input remains live for exactly length bytes; the
    // output is a nonnull, uniquely owned 32-byte array. CommonCrypto writes
    // exactly a SHA-256 digest and retains neither pointer. Empty input is valid.
    let result = unsafe { CC_SHA256(bytes.as_ptr().cast(), length, digest.as_mut_ptr()) };
    (result == digest.as_mut_ptr()).then_some(digest)
}

#[cfg(test)]
mod digest_tests {
    use super::*;
    #[test]
    fn native_digest_matches_independent_portable_backend_at_block_boundaries() {
        assert_eq!(
            sha256(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            sha256(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        for length in [0, 1, 55, 56, 63, 64, 65, 127, 128, 129, 1024, 1_048_577] {
            let input: Vec<u8> = (0..length).map(|i| (i * 73 + 17) as u8).collect();
            assert_eq!(sha256(&input), format!("{:x}", Sha256::digest(&input)));
        }
        // Exact representability boundary of the system API, without a4GiB allocation.
        assert!(u32::try_from(u32::MAX as u64).is_ok());
        assert!(u32::try_from(u32::MAX as u64 + 1).is_err());
    }
}

#[derive(Debug, Clone)]
struct Entry {
    count: u64,
    aa: String,
    line: u64,
}
#[derive(Debug)]
pub struct Sample {
    entries: BTreeMap<String, Entry>,
    pub depth: u64,
    pub rows: usize,
    pub path: PathBuf,
    pub sha256: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MergedRow {
    pub nt: String,
    pub aa: String,
    pub count1: u64,
    pub count2: u64,
}
#[derive(Debug, Serialize)]
pub struct InputSummary {
    pub path: PathBuf,
    pub sha256: String,
    pub row_count: usize,
    pub unique_clonotype_count: usize,
    pub depth: u64,
}
#[derive(Debug, Serialize)]
pub struct Provenance {
    pub command: &'static str,
    pub command_version: &'static str,
    pub inputs: [InputSummary; 2],
    pub union_clonotype_count: usize,
    pub sample_exclusive_clonotype_counts: [usize; 2],
    pub shared_clonotype_count: usize,
    pub validation_policy: &'static str,
    pub count_conservation: bool,
    pub output_path: PathBuf,
    pub output_sha256: String,
}

fn column(headers: &csv::StringRecord, aliases: &[&str], required: bool) -> Result<Option<usize>> {
    let matches: Vec<_> = headers
        .iter()
        .enumerate()
        .filter(|(_, h)| aliases.contains(h))
        .collect();
    match matches.as_slice() {
        [] if required => Err(error(format!("malformed header: missing {}", aliases[0]))),
        [] => Ok(None),
        [(index, _)] => Ok(Some(*index)),
        _ => Err(error(format!(
            "malformed header: ambiguous aliases for {}",
            aliases[0]
        ))),
    }
}

/// Parse the same bytes that are hashed, avoiding a separate hash/read race.
pub fn read_sample(path: &Path) -> Result<Sample> {
    let path = path.canonicalize().map_err(|e| context("input path", e))?;
    let bytes = fs::read(&path).map_err(|e| context("read input", e))?;
    parse_sample(&bytes, path)
}

pub fn parse_sample(bytes: &[u8], path: PathBuf) -> Result<Sample> {
    let text = std::str::from_utf8(bytes).map_err(|e| context("input must be UTF-8", e))?;
    if text.lines().any(|line| line.is_empty()) {
        return Err(error(format!(
            "{}: blank physical records are not permitted",
            path.display()
        )));
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .quoting(false)
        .from_reader(bytes);
    let headers = reader
        .headers()
        .map_err(|e| context("malformed header", e))?
        .clone();
    let mut seen = BTreeSet::new();
    if headers.iter().any(|h| h.is_empty() || !seen.insert(h)) {
        return Err(error("malformed header: empty or duplicate column label"));
    }
    let nt_col = column(&headers, &["ntCDR3", "N. Seq. CDR3"], true)?.unwrap();
    let count_col = column(&headers, &["count", "Clone count"], true)?.unwrap();
    let aa_col = column(&headers, &["AACDR3", "AA. Seq. CDR3"], false)?;
    let mut entries: BTreeMap<String, Entry> = BTreeMap::new();
    let mut depth = 0_u64;
    for record in reader.records() {
        let record =
            record.map_err(|e| context(&format!("{}: invalid TSV record", path.display()), e))?;
        let line = record.position().map_or(0, |p| p.line());
        let nt = &record[nt_col];
        if nt.trim().is_empty() {
            return Err(error(format!(
                "{} row {line}: empty nucleotide key",
                path.display()
            )));
        }
        if let Some(first) = entries.get(nt) {
            return Err(error(format!(
                "{}: duplicate nucleotide key {nt:?} at rows {} and {line}",
                path.display(),
                first.line
            )));
        }
        let raw = &record[count_col];
        if raw.is_empty() || !raw.bytes().all(|b| b.is_ascii_digit()) {
            return Err(error(format!(
                "{} row {line}: invalid count {raw:?}; require nonnegative integer ASCII digits",
                path.display()
            )));
        }
        let count = raw.parse::<u64>().map_err(|_| {
            error(format!(
                "{} row {line}: count overflow {raw:?} exceeds u64",
                path.display()
            ))
        })?;
        depth = depth.checked_add(count).ok_or_else(|| {
            error(format!(
                "{} row {line}: sample depth overflow",
                path.display()
            ))
        })?;
        entries.insert(
            nt.to_owned(),
            Entry {
                count,
                aa: aa_col.map_or("", |i| &record[i]).to_owned(),
                line,
            },
        );
    }
    if depth == 0 {
        return Err(error(format!(
            "{}: sample total depth is zero",
            path.display()
        )));
    }
    Ok(Sample {
        rows: entries.len(),
        entries,
        depth,
        path,
        sha256: sha256(bytes),
    })
}

pub fn merge(first: &Sample, second: &Sample) -> Result<Vec<MergedRow>> {
    let mut left = first.entries.iter().peekable();
    let mut right = second.entries.iter().peekable();
    let mut rows = Vec::with_capacity(first.entries.len().max(second.entries.len()));
    loop {
        let ordering = match (left.peek(), right.peek()) {
            (Some((a, _)), Some((b, _))) => a.cmp(b),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => break,
        };
        let (key, a, b) = match ordering {
            std::cmp::Ordering::Less => {
                let (key, value) = left.next().unwrap();
                (key, Some(value), None)
            }
            std::cmp::Ordering::Greater => {
                let (key, value) = right.next().unwrap();
                (key, None, Some(value))
            }
            std::cmp::Ordering::Equal => {
                let (key, a) = left.next().unwrap();
                let (_, b) = right.next().unwrap();
                (key, Some(a), Some(b))
            }
        };
        if let (Some(a), Some(b)) = (a, b) {
            if !a.aa.is_empty() && !b.aa.is_empty() && a.aa != b.aa {
                return Err(error(format!("conflicting amino-acid annotations for nucleotide key {key:?}: {} row {} {:?} versus {} row {} {:?}", first.path.display(), a.line, a.aa, second.path.display(), b.line, b.aa)));
            }
        }
        let aa = a
            .filter(|e| !e.aa.is_empty())
            .or(b)
            .map_or("", |e| &e.aa)
            .to_owned();
        rows.push(MergedRow {
            nt: key.clone(),
            aa,
            count1: a.map_or(0, |e| e.count),
            count2: b.map_or(0, |e| e.count),
        });
    }
    check_conservation(&rows, [first.depth, second.depth])?;
    Ok(rows)
}

pub fn check_conservation(rows: &[MergedRow], depths: [u64; 2]) -> Result<()> {
    let totals = rows.iter().try_fold([0_u64; 2], |acc, row| {
        Some([
            acc[0].checked_add(row.count1)?,
            acc[1].checked_add(row.count2)?,
        ])
    });
    if totals != Some(depths) {
        return Err(error(
            "count conservation failed: merged depths differ or overflow",
        ));
    }
    Ok(())
}

fn encode(rows: &[MergedRow]) -> Result<Vec<u8>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Vec::new());
    writer
        .write_record(["ntCDR3", "AACDR3", "count_1", "count_2"])
        .map_err(|e| context("encode header", e))?;
    for row in rows {
        writer
            .write_record([
                &row.nt,
                &row.aa,
                &row.count1.to_string(),
                &row.count2.to_string(),
            ])
            .map_err(|e| context("encode row", e))?;
    }
    writer.into_inner().map_err(|e| context("finish TSV", e))
}

fn write_synced(path: &Path, bytes: &[u8]) -> Result<()> {
    let mut file = File::create(path).map_err(|e| context("create staged output", e))?;
    file.write_all(bytes)
        .and_then(|()| file.sync_all())
        .map_err(|e| context("write staged output", e))
}

/// Publish a pair of files with one atomic, no-replace directory rename.
/// The parent must exist. No overwrite fallback is used on unsupported platforms.
pub fn validate(first: &Path, second: &Path, output_dir: &Path) -> Result<Provenance> {
    validate_inner(first, second, output_dir, |_, _, _| Ok(()))
}

fn validate_inner(
    first: &Path,
    second: &Path,
    output_dir: &Path,
    hook: impl Fn(u8, &Path, &Path) -> Result<()>,
) -> Result<Provenance> {
    let a = read_sample(first)?;
    let b = read_sample(second)?;
    let rows = merge(&a, &b)?;
    let parent = output_dir
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent = parent
        .canonicalize()
        .map_err(|e| context("output parent must exist", e))?;
    let name = output_dir
        .file_name()
        .filter(|n| !n.is_empty())
        .ok_or_else(|| error("output directory needs a new final name"))?;
    let output = parent.join(name);
    if fs::symlink_metadata(&output).is_ok() {
        return Err(error(format!(
            "refusing existing output: {}",
            output.display()
        )));
    }
    let merged = encode(&rows)?;
    let shared = a
        .entries
        .keys()
        .filter(|k| b.entries.contains_key(*k))
        .count();
    let input_summary = |s: &Sample| InputSummary {
        path: s.path.clone(),
        sha256: s.sha256.clone(),
        row_count: s.rows,
        unique_clonotype_count: s.entries.len(),
        depth: s.depth,
    };
    let provenance = Provenance {
        command: "clonoshift validate",
        command_version: VERSION,
        inputs: [input_summary(&a), input_summary(&b)],
        union_clonotype_count: rows.len(),
        sample_exclusive_clonotype_counts: [a.entries.len() - shared, b.entries.len() - shared],
        shared_clonotype_count: shared,
        validation_policy: POLICY,
        count_conservation: true,
        output_path: output.join("merged.tsv"),
        output_sha256: sha256(&merged),
    };
    let mut summary =
        serde_json::to_vec_pretty(&provenance).map_err(|e| context("encode provenance", e))?;
    summary.push(b'\n');
    let staging = tempfile::Builder::new()
        .prefix(".clonoshift-stage-")
        .tempdir_in(&parent)
        .map_err(|e| context("create staging directory", e))?;
    write_synced(&staging.path().join("merged.tsv"), &merged)?;
    hook(1, staging.path(), &output)?;
    write_synced(&staging.path().join("provenance.json"), &summary)?;
    hook(2, staging.path(), &output)?;
    File::open(staging.path())
        .and_then(|f| f.sync_all())
        .map_err(|e| context("sync staging directory", e))?;
    publish(staging.path(), &output)
        .map_err(|e| context("atomic output publication (no replacement)", e))?;
    // TempDir's old path no longer exists; its cleanup cannot remove the published destination.
    Ok(provenance)
}

#[cfg(any(target_os = "macos", target_os = "linux"))]
pub fn publish(source: &Path, destination: &Path) -> std::io::Result<()> {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;
    let src = CString::new(source.as_os_str().as_bytes())?;
    let dst = CString::new(destination.as_os_str().as_bytes())?;
    // SAFETY: CString pointers are NUL-terminated and live for the duration of the call.
    #[cfg(target_os = "macos")]
    let result = unsafe {
        libc::renameatx_np(
            libc::AT_FDCWD,
            src.as_ptr(),
            libc::AT_FDCWD,
            dst.as_ptr(),
            libc::RENAME_EXCL,
        )
    };
    #[cfg(target_os = "linux")]
    let result = unsafe {
        libc::renameat2(
            libc::AT_FDCWD,
            src.as_ptr(),
            libc::AT_FDCWD,
            dst.as_ptr(),
            libc::RENAME_NOREPLACE,
        )
    };
    if result == 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error())
    }
}
#[cfg(not(any(target_os = "macos", target_os = "linux")))]
pub fn publish(_: &Path, _: &Path) -> std::io::Result<()> {
    Err(std::io::Error::new(
        std::io::ErrorKind::Unsupported,
        "atomic no-replace publication supports macOS and Linux only",
    ))
}

#[cfg(test)]
mod atomic_tests {
    use super::*;
    #[test]
    fn no_replace_race_and_failure_preserve_staging() {
        let root = tempfile::tempdir().unwrap();
        let stage = root.path().join("stage");
        let destination = root.path().join("out");
        fs::create_dir(&stage).unwrap();
        fs::write(stage.join("merged.tsv"), b"staged").unwrap();
        fs::create_dir(&destination).unwrap();
        assert!(publish(&stage, &destination).is_err());
        assert!(destination.read_dir().unwrap().next().is_none());
        assert_eq!(fs::read(stage.join("merged.tsv")).unwrap(), b"staged");
        fs::remove_dir(&destination).unwrap();
        publish(&stage, &destination).unwrap();
        assert!(!stage.exists());
        assert_eq!(fs::read(destination.join("merged.tsv")).unwrap(), b"staged");
    }

    #[test]
    fn transaction_failures_clean_staging_without_partial_publication() {
        let root = tempfile::tempdir().unwrap();
        let first = root.path().join("a.tsv");
        let second = root.path().join("b.tsv");
        fs::write(&first, b"ntCDR3\tcount\nA\t1\n").unwrap();
        fs::write(&second, b"ntCDR3\tcount\nB\t2\n").unwrap();
        for phase in [1, 2] {
            let out = root.path().join(format!("out{phase}"));
            let result = validate_inner(&first, &second, &out, |p, stage, output| {
                assert!(!output.exists());
                assert!(stage.join("merged.tsv").is_file());
                if p == phase {
                    Err(error("injected prepublication failure"))
                } else {
                    Ok(())
                }
            });
            assert!(result.unwrap_err().to_string().contains("injected"));
            assert!(!out.exists());
            assert!(!root.path().read_dir().unwrap().any(|p| p
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".clonoshift-stage-")));
        }
        let out = root.path().join("raced");
        let result = validate_inner(&first, &second, &out, |phase, _, destination| {
            if phase == 2 {
                fs::create_dir(destination).unwrap();
            }
            Ok(())
        });
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("atomic output publication"));
        assert!(out.read_dir().unwrap().next().is_none());
        assert!(!root.path().read_dir().unwrap().any(|p| p
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".clonoshift-stage-")));
    }
}
