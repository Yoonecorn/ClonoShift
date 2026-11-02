// SPDX-License-Identifier: GPL-3.0-or-later
use clonoshift::acquisition::Acquisition;
use clonoshift::autofit::{self, AutoFitConfig};
use clonoshift::fitting::evaluate_fixed;
use clonoshift::observation::{ObservedCountPair, PairTable};
fn table(xs: &[[u64; 3]]) -> PairTable {
    PairTable::new(
        &xs.iter()
            .map(|r| ObservedCountPair {
                count1: r[0],
                count2: r[1],
                multiplicity: r[2],
            })
            .collect::<Vec<_>>(),
    )
    .unwrap()
}
#[test]
fn poisson_total_is_sufficient_for_c2_and_parameter_likelihood_differences() {
    let pairs = table(&[
        [1, 0, 8],
        [0, 1, 5],
        [2, 1, 4],
        [1, 2, 2],
        [0, 9, 1],
        [4, 5, 3],
    ]);
    let totals = table(&[[1, 0, 13], [3, 0, 6], [9, 0, 4]]);
    let mut constants = vec![];
    for (a, x) in [(-2., -8.), (-1.7, -6.), (-2.4, -5.)] {
        let original =
            evaluate_fixed(&pairs, [100., 200.], Acquisition::Poisson, a, x, 257).unwrap();
        let sufficient =
            evaluate_fixed(&totals, [300., 0.], Acquisition::Poisson, a, x, 257).unwrap();
        assert!((original.objective.c2 - sufficient.objective.c2).abs() < 2e-13);
        assert!((original.objective.c1 - sufficient.objective.c1).abs() < 2e-13);
        constants.push(
            original.objective.negative_log_likelihood
                - sufficient.objective.negative_log_likelihood,
        );
    }
    assert!((constants[0] - constants[1]).abs() < 2e-12);
    assert!((constants[0] - constants[2]).abs() < 2e-12);
}
#[test]
fn refuses_empty_unrepresentable_and_out_of_count_domain_explicitly() {
    let empty = table(&[]);
    assert_eq!(
        autofit::fit(&empty, [10., 10.], &AutoFitConfig::default()).outcome,
        "Refused"
    );
    let high = table(&[[1_000_001, 0, 1]]);
    let r = autofit::fit(&high, [100., 100.], &AutoFitConfig::default());
    assert!(r.reason.contains("exceeds"));
    assert!(r.candidate.is_none());
    let high = table(&[[u64::MAX, 1, 1]]);
    assert!(autofit::fit(&high, [100., 100.], &AutoFitConfig::default())
        .reason
        .contains("overflow"));
    let ordinary = table(&[[1, 0, 3]]);
    for depths in [[0., 10.], [f64::INFINITY, 1.], [f64::NAN, 1.], [1e16, 1.]] {
        let r = autofit::fit(&ordinary, depths, &AutoFitConfig::default());
        assert_eq!(r.outcome, "Refused");
        assert!(r.candidate.is_none());
    }
}
#[test]
fn multistart_report_is_deterministic_and_swaps_preserve_fit() {
    let original = table(&[
        [1, 0, 8],
        [0, 1, 7],
        [1, 1, 6],
        [2, 1, 5],
        [1, 2, 4],
        [3, 3, 2],
        [5, 6, 2],
    ]);
    let swapped = table(&[
        [0, 1, 8],
        [1, 0, 7],
        [1, 1, 6],
        [1, 2, 5],
        [2, 1, 4],
        [3, 3, 2],
        [6, 5, 2],
    ]);
    let c = AutoFitConfig {
        log_minimum_bounds: Some([-12., -0.01]),
        ..AutoFitConfig::default()
    };
    let a = autofit::fit(&original, [100., 100.], &c);
    assert_eq!(a, autofit::fit(&original, [100., 100.], &c));
    let b = autofit::fit(&swapped, [100., 100.], &c);
    assert_eq!(a.outcome, b.outcome);
    assert_eq!(a.grids, b.grids);
    assert_eq!(
        a.diagnostic_candidate.as_ref().map(|x| x.exponent),
        b.diagnostic_candidate.as_ref().map(|x| x.exponent)
    );
    assert_eq!(
        a.diagnostic_candidate.as_ref().map(|x| x.minimum),
        b.diagnostic_candidate.as_ref().map(|x| x.minimum)
    );
    assert!(a.parameter_uncertainty.contains("Unavailable"));
}
#[test]
fn singletons_are_explicitly_uninformative_for_allocation_dispersion() {
    let singleton = table(&[[1, 0, 100], [0, 1, 100]]);
    let c = AutoFitConfig {
        grid_points: 65,
        starts: vec![-2.],
        ..AutoFitConfig::default()
    };
    let r = autofit::fit(&singleton, [100., 100.], &c);
    let s = r.split_pearson.unwrap();
    assert_eq!(s.informative_observations, 0);
    assert_eq!(s.pearson_mean, None);
    assert_eq!(r.outcome, "Refused");
    assert!(r.candidate.is_none());
}
