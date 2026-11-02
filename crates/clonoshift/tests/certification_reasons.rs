// SPDX-License-Identifier: GPL-3.0-or-later
use clonoshift::acquisition::Acquisition;
use clonoshift::certification::CertificationReason as Reason;
use clonoshift::observation::SampleConfiguration;
use clonoshift::response::{JointPrior, JointState, ResponseError, ResponseModel};

fn model(states: &[[f64; 3]], depths: [f64; 2], symmetric: bool) -> ResponseModel {
    let prior = JointPrior::new(
        &states
            .iter()
            .map(|s| JointState {
                frequency1: s[0],
                frequency2: s[1],
                mass: s[2],
            })
            .collect::<Vec<_>>(),
        symmetric,
    )
    .unwrap();
    let samples = depths.map(|depth| SampleConfiguration {
        acquisition: Acquisition::Poisson,
        depth,
    });
    ResponseModel::new(prior, samples[0], samples[1], 0.3, 1.0).unwrap()
}
fn check(m: &ResponseModel, counts: [u64; 2], expected: Reason) {
    let failure = m
        .infer_certified_detailed(counts[0], counts[1])
        .unwrap_err();
    assert_eq!(failure.error, ResponseError::PrecisionLoss);
    assert_eq!(failure.reason, Some(expected));
    assert_eq!(
        m.infer_certified(counts[0], counts[1]).unwrap_err(),
        ResponseError::PrecisionLoss
    );
    assert_eq!(
        failure,
        m.infer_certified_detailed(counts[0], counts[1])
            .unwrap_err()
    );
}
#[test]
fn anchor_failure_reports_reason() {
    check(
        &model(
            &[[0.9, 0.9, 0.5], [1.0, 1.0, 0.5]],
            [200000.0, 200000.0],
            true,
        ),
        [1, 1],
        Reason::RetainedEvidenceAnchorUnavailable,
    );
}
#[test]
fn loose_cumulative_bound_reports_reason() {
    check(
        &model(
            &[[0.70365, 0.70365, 0.5], [0.7037, 0.7037, 0.5]],
            [200000.0, 200000.0],
            true,
        ),
        [1, 0],
        Reason::CumulativeOmissionBoundTooLoose,
    );
}
#[test]
fn ambiguous_interval_reports_reason() {
    check(
        &model(
            &[[1e-3, 1e-5, 0.0125], [1e-5, 1e-5, 0.4875], [1.0, 1.0, 0.5]],
            [0.0, 200000.0],
            false,
        ),
        [0, 1],
        Reason::AmbiguousIntervalCertification,
    );
}
#[test]
fn unsupported_certification_configuration_reports_reason() {
    check(
        &model(
            &[[1e-5, 1e-5, 0.5], [1.0, 1.0, 0.5]],
            [200000.0, 200000.0],
            true,
        ),
        [1_000_001, 1],
        Reason::InvalidCertificationConfiguration,
    );
}

#[test]
fn all_reason_codes_have_stable_machine_readable_spelling() {
    use clonoshift::analysis::{error_json, AnalysisError, State};
    for (reason, expected) in [
        (
            Reason::RetainedEvidenceAnchorUnavailable,
            "RetainedEvidenceAnchorUnavailable",
        ),
        (
            Reason::CumulativeOmissionBoundTooLoose,
            "CumulativeOmissionBoundTooLoose",
        ),
        (
            Reason::UnsupportedArithmeticDomain,
            "UnsupportedArithmeticDomain",
        ),
        (
            Reason::AmbiguousIntervalCertification,
            "AmbiguousIntervalCertification",
        ),
        (
            Reason::InvalidCertificationConfiguration,
            "InvalidCertificationConfiguration",
        ),
        (
            Reason::InternalCertificationFailure,
            "InternalCertificationFailure",
        ),
    ] {
        let error = AnalysisError {
            state: State::UnresolvedResult,
            message: reason.message().into(),
            certification_reason: Some(reason),
        };
        let json: serde_json::Value = serde_json::from_str(&error_json(&error)).unwrap();
        assert_eq!(json["classification"], "PrecisionLoss");
        assert_eq!(json["reason_code"], expected);
    }
}
