// SPDX-License-Identifier: GPL-3.0-or-later
//! Pinned NoisET compatibility only. Retains first-category conditioning and
//! defective unseen-mass C2. Never used by corrected research modules.
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const REVISION: &str = "2b6a45fcd26e7f9eefb3004d81725355cf6d66c9";
pub const POINTS: usize = 1200;
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pair {
    pub counts: [u64; 2],
    pub multiplicity: u64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Data {
    pub pairs: Vec<Pair>,
    pub depths: [f64; 2],
    pub observations: u64,
}
impl Data {
    pub fn new(counts: &[[u64; 2]]) -> Result<Self, String> {
        let mut grouped = BTreeMap::new();
        let mut depths = [0_u64; 2];
        if counts.is_empty() {
            return Err("empty legacy observations".into());
        }
        for &pair in counts {
            if pair == [0, 0] || pair.iter().any(|&n| n > 1_000_000) {
                return Err(
                    "legacy supported observations require nonzero pair and counts <= 1000000"
                        .into(),
                );
            }
            for i in 0..2 {
                depths[i] = depths[i].checked_add(pair[i]).ok_or("depth overflow")?;
            }
            *grouped.entry(pair).or_insert(0) += 1;
        }
        if depths.contains(&0) || depths.iter().any(|&d| d > (1_u64 << 53)) {
            return Err(
                "legacy supported depths must be positive and exactly representable".into(),
            );
        }
        Ok(Self {
            pairs: grouped
                .into_iter()
                .map(|(counts, multiplicity)| Pair {
                    counts,
                    multiplicity,
                })
                .collect(),
            depths: depths.map(|d| d as f64),
            observations: counts.len() as u64,
        })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Evaluation {
    pub parameters: [f64; 2],
    pub log_frequency: Vec<f64>,
    pub log_prior: Vec<f64>,
    pub log_normalizer: f64,
    pub first_category: [u64; 2],
    pub first_probability: f64,
    pub explicit_zero_probability: f64,
    pub pair_probabilities: Vec<f64>,
    pub conditioned_probabilities: Vec<f64>,
    pub mean_nll: f64,
    pub c1: f64,
    pub c2: f64,
    pub diversity: u64,
}
fn integral(x: &[f64], y: &[f64]) -> f64 {
    x.windows(2)
        .zip(y.windows(2))
        .map(|(x, y)| (x[1] - x[0]) * 0.5 * (y[0] + y[1]))
        .sum()
}
pub fn evaluate(data: &Data, parameters: [f64; 2]) -> Result<Evaluation, String> {
    let [a, b] = parameters;
    if !a.is_finite() || !b.is_finite() || b >= 0.0 || 10_f64.powf(b) == 0.0 {
        return Err("unrepresentable legacy prior support".into());
    }
    // This public struct is serializable; revalidate rather than trusting fields.
    if data.pairs.is_empty()
        || data.observations == 0
        || data.depths.iter().any(|d| !d.is_finite() || *d <= 0.0)
        || data.pairs.iter().any(|p| {
            p.multiplicity == 0 || p.counts == [0, 0] || p.counts.iter().any(|&n| n > 1_000_000)
        })
        || data
            .pairs
            .iter()
            .map(|p| p.multiplicity as u128)
            .sum::<u128>()
            != data.observations as u128
    {
        return Err("invalid legacy data".into());
    }
    if data.pairs.windows(2).any(|p| p[0].counts >= p[1].counts) {
        return Err("legacy pairs must be uniquely sorted".into());
    }
    for axis in 0..2 {
        let sum = data
            .pairs
            .iter()
            .try_fold(0_u128, |s, p| {
                s.checked_add(p.counts[axis] as u128 * p.multiplicity as u128)
            })
            .ok_or("weighted depth overflow")?;
        if sum > (1_u128 << 53) || data.depths[axis] != sum as f64 {
            return Err("legacy weighted depth mismatch".into());
        }
    }
    evaluate_inner(data, parameters)
}
/// Fixed-depth component diagnostic for compatibility checks. The file workflow
/// always uses count-conserved depths; these caller depths are never fitted.
pub fn evaluate_at_depths(
    data: &Data,
    parameters: [f64; 2],
    depths: [f64; 2],
) -> Result<Evaluation, String> {
    evaluate(data, parameters)?;
    if depths.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err("invalid diagnostic depths".into());
    }
    let mut data = data.clone();
    data.depths = depths;
    evaluate_inner(&data, parameters)
}
fn evaluate_inner(data: &Data, parameters: [f64; 2]) -> Result<Evaluation, String> {
    let [a, b] = parameters;
    let log_frequency: Vec<_> = (0..POINTS)
        .map(|i| {
            10_f64
                .powf(b + (0.0 - b) * i as f64 / (POINTS - 1) as f64)
                .ln()
        })
        .collect();
    let z = integral(
        &log_frequency,
        &log_frequency
            .iter()
            .map(|&f| ((a + 1.0) * f).exp())
            .collect::<Vec<_>>(),
    )
    .ln();
    let log_prior: Vec<_> = log_frequency.iter().map(|&f| a * f - z).collect();
    let max = data.pairs.iter().flat_map(|p| p.counts).max().unwrap() as usize;
    let mut factorial = vec![0.0; max + 1];
    for i in 1..=max {
        factorial[i] = factorial[i - 1] + (i as f64).ln();
    }
    let first_category = [0, 1].map(|i| data.pairs.iter().map(|p| p.counts[i]).min().unwrap());
    // Reuse identical transcendental results; keep every arithmetic expression
    // and trapezoid accumulation in its original order. In particular exp().ln()
    // is intentional: pinned compatibility includes upstream PMF underflow.
    let frequencies: Vec<_> = log_frequency.iter().map(|&f| f.exp()).collect();
    let means: [Vec<f64>; 2] =
        std::array::from_fn(|axis| frequencies.iter().map(|&f| data.depths[axis] * f).collect());
    let log_means: [Vec<f64>; 2] =
        std::array::from_fn(|axis| means[axis].iter().map(|&m| m.ln()).collect());
    let mut pmfs: [BTreeMap<u64, Vec<f64>>; 2] = std::array::from_fn(|_| BTreeMap::new());
    for axis in 0..2 {
        // Explicit zero is required even when no observed category is zero.
        for count in data
            .pairs
            .iter()
            .map(|p| p.counts[axis])
            .chain(std::iter::once(0))
        {
            pmfs[axis].entry(count).or_insert_with(|| {
                means[axis]
                    .iter()
                    .zip(&log_means[axis])
                    .map(|(&m, &log_m)| {
                        (count as f64 * log_m - m - factorial[count as usize])
                            .exp()
                            .ln()
                    })
                    .collect()
            });
        }
    }
    let pair_integrals = |counts: [u64; 2]| {
        let l1 = &pmfs[0][&counts[0]];
        let l2 = &pmfs[1][&counts[1]];
        let y: Vec<_> = log_frequency
            .iter()
            .zip(&log_prior)
            .zip(l1.iter().zip(l2))
            .map(|((&f, &rho), (&l1, &l2))| (rho + f + l1 + l2).exp())
            .collect();
        let mass = integral(&log_frequency, &y);
        let moment = integral(
            &log_frequency,
            &y.iter()
                .zip(&frequencies)
                .map(|(&p, &f)| p * f)
                .collect::<Vec<_>>(),
        );
        (mass, moment)
    };
    let (p0, m0) = pair_integrals(first_category);
    let explicit_zero_probability = pair_integrals([0, 0]).0;
    let mut pair_probabilities = Vec::new();
    let mut conditioned_probabilities = Vec::new();
    let mut nll = 0.0;
    let mut observed_mass = 0.0;
    for p in &data.pairs {
        let (prob, moment) = pair_integrals(p.counts);
        let conditioned = prob / (1.0 - p0);
        // Upstream silently assigns log-loss zero to underflowed probabilities.
        if conditioned > 0.0 {
            nll -= p.multiplicity as f64 * conditioned.ln();
        }
        if prob > 0.0 {
            observed_mass += p.multiplicity as f64 * moment / prob;
        }
        pair_probabilities.push(prob);
        conditioned_probabilities.push(conditioned);
    }
    let inferred = data.observations as f64 / (1.0 - p0);
    let mean = integral(
        &log_frequency,
        &log_frequency
            .iter()
            .zip(&log_prior)
            .map(|(&f, &rho)| (rho + 2.0 * f).exp())
            .collect::<Vec<_>>(),
    );
    let c1 = mean.ln() - ((1.0 - p0) / (data.observations as f64)).ln();
    let c2 = (inferred * p0 * m0 + observed_mass).ln();
    if log_prior.iter().any(|v| !v.is_finite())
        || !z.is_finite()
        || !(0.0..1.0).contains(&p0)
        || !nll.is_finite()
        || !c1.is_finite()
        || !c2.is_finite()
        || !inferred.is_finite()
        || inferred >= u64::MAX as f64
    {
        return Err("legacy arithmetic outside supported finite range".into());
    }
    Ok(Evaluation {
        parameters,
        log_frequency,
        log_prior,
        log_normalizer: z,
        first_category,
        first_probability: p0,
        explicit_zero_probability,
        pair_probabilities,
        conditioned_probabilities,
        mean_nll: nll / data.observations as f64,
        c1,
        c2,
        diversity: inferred as u64,
    })
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Iteration {
    pub parameters: [f64; 2],
    pub mean_nll: f64,
    pub c2: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Fit {
    pub success: bool,
    pub status: String,
    pub iterations: usize,
    pub initial: [f64; 2],
    pub trace: Vec<Iteration>,
    pub evaluation: Evaluation,
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
/// Equality SQP with a BFGS Lagrangian Hessian. Same legacy objective,
/// initialization, no parameter bounds, equality, ftol and iteration limit.
/// Algorithm/trace differs from SciPy SLSQP and is separately qualified.
pub fn fit(data: &Data, initial: [f64; 2], max_iterations: usize) -> Result<Fit, String> {
    let mut x = initial;
    let mut e = evaluate(data, x)?;
    let mut h: [[f64; 2]; 2] = [[1.0, 0.0], [0.0, 1.0]];
    let mut trace = Vec::new();
    let mut success = false;
    let mut status = "iteration limit";
    let gradient = |x: [f64; 2], e: &Evaluation| -> Result<([f64; 2], [f64; 2]), String> {
        let mut g = [0.0; 2];
        let mut c = [0.0; 2];
        for i in 0..2 {
            let step = 1e-5;
            let mut l = x;
            let mut r = x;
            l[i] -= step;
            r[i] += step;
            let le = evaluate(data, l)?;
            let re = evaluate(data, r)?;
            g[i] = (re.mean_nll - le.mean_nll) / (2.0 * step);
            c[i] = (re.c2 - le.c2) / (2.0 * step);
        }
        let _ = e;
        Ok((g, c))
    };
    for _ in 0..max_iterations.min(200) {
        let (g, c) = gradient(x, &e)?;
        let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
        if !det.is_finite() || det <= 1e-15 {
            h = [[1.0, 0.0], [0.0, 1.0]];
        }
        let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
        let inverse = |v: [f64; 2]| {
            [
                (h[1][1] * v[0] - h[0][1] * v[1]) / det,
                (-h[1][0] * v[0] + h[0][0] * v[1]) / det,
            ]
        };
        let hg = inverse(g);
        let hc = inverse(c);
        let denom = dot(c, hc);
        if !denom.is_finite() || denom <= 1e-20 {
            status = "singular equality";
            break;
        }
        let lambda = (dot(c, hg) - e.c2) / denom;
        let p = [-hg[0] + lambda * hc[0], -hg[1] + lambda * hc[1]];
        let stationarity = dot(
            [g[0] - lambda * c[0], g[1] - lambda * c[1]],
            [g[0] - lambda * c[0], g[1] - lambda * c[1]],
        )
        .sqrt();
        if e.c2.abs() < 1e-6 && stationarity < 1e-7 {
            success = true;
            status = "converged";
            break;
        }
        let penalty = 10.0_f64.max(2.0 * lambda.abs());
        let merit = e.mean_nll + penalty * e.c2.abs();
        let mut step = 1.0;
        let mut trial = None;
        for _ in 0..40 {
            let nx = [x[0] + step * p[0], x[1] + step * p[1]];
            if let Ok(ne) = evaluate(data, nx) {
                if ne.mean_nll + penalty * ne.c2.abs() < merit {
                    trial = Some((nx, ne));
                    break;
                }
            }
            step *= 0.5;
        }
        let Some((nx, ne)) = trial else {
            status = "line search failed";
            break;
        };
        let (ng, nc) = gradient(nx, &ne)?;
        let s = [nx[0] - x[0], nx[1] - x[1]];
        let y = [
            ng[0] - g[0] - lambda * (nc[0] - c[0]),
            ng[1] - g[1] - lambda * (nc[1] - c[1]),
        ];
        let hs = [dot(h[0], s), dot(h[1], s)];
        let sy = dot(s, y);
        let shs = dot(s, hs);
        if sy > 1e-12 && shs > 1e-12 {
            for i in 0..2 {
                for j in 0..2 {
                    h[i][j] += y[i] * y[j] / sy - hs[i] * hs[j] / shs;
                }
            }
        }
        x = nx;
        e = ne;
        trace.push(Iteration {
            parameters: x,
            mean_nll: e.mean_nll,
            c2: e.c2,
        });
    }
    Ok(Fit {
        success,
        status: status.into(),
        iterations: trace.len(),
        initial,
        trace,
        evaluation: e,
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Posterior {
    pub posterior: Vec<f64>,
    pub tail: f64,
    pub mean: f64,
    pub mode: f64,
    pub lower: Option<f64>,
    pub median: Option<f64>,
    pub upper: Option<f64>,
}
pub fn response_prior(alpha: f64, scale: f64, step: f64) -> Result<Vec<f64>, String> {
    if !alpha.is_finite()
        || !(0.0..=1.0).contains(&alpha)
        || !scale.is_finite()
        || scale <= 0.0
        || !step.is_finite()
        || step <= 0.0
    {
        return Err("invalid legacy response prior".into());
    }
    let lambda = -step / scale;
    let z = 2.0 * ((251.0 * lambda).exp() - 1.0) / (lambda.exp() - 1.0) - 1.0;
    let mut p: Vec<_> = (0..501)
        .map(|i| alpha * (lambda * (i as f64 - 250.0).abs()).exp() / z)
        .collect();
    p[250] += 1.0 - alpha;
    if !z.is_finite() || z <= 0.0 || p.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err("unrepresentable legacy response prior".into());
    }
    Ok(p)
}
pub fn summarize(shifts: &[f64], prior: &[f64], likelihood: &[f64]) -> Result<Posterior, String> {
    if shifts.is_empty()
        || shifts.len() != prior.len()
        || shifts.len() != likelihood.len()
        || shifts.iter().any(|v| !v.is_finite())
        || shifts.windows(2).any(|w| w[1] <= w[0])
        || prior
            .iter()
            .chain(likelihood)
            .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err("invalid legacy posterior arrays".into());
    }
    let mut posterior: Vec<_> = prior.iter().zip(likelihood).map(|(p, l)| p * l).collect();
    let total: f64 = posterior.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err("zero legacy posterior evidence".into());
    }
    posterior.iter_mut().for_each(|p| *p /= total);
    let mut cdf = Vec::new();
    let mut sum = 0.0;
    for p in &posterior {
        sum += p;
        cdf.push(sum);
    }
    let mut back = vec![0.0; posterior.len()];
    sum = 0.0;
    for i in (0..posterior.len()).rev() {
        sum += posterior[i];
        back[i] = sum;
    }
    let zero = shifts
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        .unwrap()
        .0;
    let mode = posterior
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1).then_with(|| b.0.cmp(&a.0)))
        .unwrap()
        .0;
    let quant = |q: f64| {
        cdf.windows(2)
            .position(|v| v[0] < q && v[1] >= q)
            .map(|i| shifts[i + 1])
    };
    let middle: Vec<_> = (0..shifts.len())
        .filter(|&i| cdf[i] >= 0.5 && back[i] >= 0.5)
        .map(|i| shifts[i])
        .collect();
    Ok(Posterior {
        tail: cdf[zero],
        mean: dot_dynamic(shifts, &posterior),
        mode: shifts[mode],
        lower: quant(0.025),
        median: if middle.is_empty() {
            None
        } else {
            Some(middle.iter().sum::<f64>() / middle.len() as f64)
        },
        upper: quant(0.975),
        posterior,
    })
}
fn dot_dynamic(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Response {
    pub parameters1: [f64; 2],
    pub parameters2: [f64; 2],
    pub shifts: Vec<f64>,
    pub pair_likelihoods: Vec<Vec<f64>>,
    pub first_likelihood: Vec<f64>,
    pub surface: Vec<Vec<f64>>,
    pub best_index: [usize; 2],
    pub alpha: f64,
    pub scale: f64,
    pub reporting_prior: Vec<f64>,
    pub summaries: Vec<Posterior>,
}
/// Aliased public NoisET Poisson response path. Preserves the extended f2>1
/// domain and the learner/reporting prior-step mismatch. No contraction test.
pub fn response(
    data: &Data,
    parameters1: [f64; 2],
    parameters2: [f64; 2],
) -> Result<Response, String> {
    let base = evaluate(data, parameters1)?;
    evaluate(data, parameters2)?;
    let x = &base.log_frequency;
    let fstep = x[1] - x[0];
    let stride = (0.1 / fstep).round_ties_even() as usize;
    if stride == 0 || stride > 10_000 {
        return Err("unsupported legacy response shift stride".into());
    }
    let step = stride as f64 * fstep;
    let shifts: Vec<_> = (0..501).map(|i| (i as f64 - 250.0) * step).collect();
    let wide_len = POINTS + 500 * stride;
    let low = x[0] - stride as f64 * 250.0 * fstep;
    let high = x[POINTS - 1] + stride as f64 * 250.0 * fstep;
    let wide: Vec<_> = (0..wide_len)
        .map(|i| low + (high - low) * i as f64 / (wide_len - 1) as f64)
        .collect();
    let max = data.pairs.iter().flat_map(|p| p.counts).max().unwrap() as usize;
    let mut fact = vec![0.0; max + 1];
    for i in 1..=max {
        fact[i] = fact[i - 1] + (i as f64).ln();
    }
    let logpmf = |n: u64, d: f64, f: f64| {
        let mean = d * f.exp();
        (n as f64 * mean.ln() - mean - fact[n as usize]).exp().ln()
    };
    let mut cats: [BTreeMap<u64, Vec<f64>>; 2] = [BTreeMap::new(), BTreeMap::new()];
    for (i, cat) in cats.iter_mut().enumerate() {
        for p in &data.pairs {
            cat.entry(p.counts[i]).or_insert_with(|| {
                (if i == 0 { x } else { &wide })
                    .iter()
                    .map(|&f| logpmf(p.counts[i], data.depths[i], f))
                    .collect()
            });
        }
    }
    let curve = |counts: [u64; 2]| -> Vec<f64> {
        let l = &cats[0][&counts[0]];
        let r = &cats[1][&counts[1]];
        (0..501)
            .map(|s| {
                let y: Vec<_> = (0..POINTS)
                    .map(|i| (base.log_prior[i] + l[i] + r[s * stride + i] + x[i]).exp())
                    .collect();
                integral(x, &y)
            })
            .collect()
    };
    // Experimental bounded parallelism: only independent observed-pair
    // curves run concurrently. Each curve and all subsequent reductions keep
    // their original arithmetic order. The default uses at most four available workers.
    let workers = match std::env::var("CLONOSHIFT_WORKERS").as_deref() {
        Ok("1") => 1,
        Ok("2") => 2,
        Ok("4") => 4,
        _ => std::thread::available_parallelism()
            .map(|n| n.get().min(4))
            .unwrap_or(1),
    };
    let pair_likelihoods: Vec<_> = if workers == 1 || data.pairs.len() < 64 {
        data.pairs.iter().map(|p| curve(p.counts)).collect()
    } else {
        let chunk_size = data.pairs.len().div_ceil(workers);
        std::thread::scope(|scope| {
            let mut chunks = data.pairs.chunks(chunk_size);
            let first = chunks.next().expect("validated nonempty pairs");
            let mut pending = Vec::new();
            let mut completed = Vec::new();
            // The calling thread owns chunk0, so requested workers includes
            // the caller; at most three additional threads can be active.
            for (offset, chunk) in chunks.enumerate() {
                let index = offset + 1;
                let curve_ref = &curve;
                let task = move || {
                    chunk
                        .iter()
                        .map(|p| curve_ref(p.counts))
                        .collect::<Vec<_>>()
                };
                match std::thread::Builder::new()
                    .name(format!("clonoshift-response-{index}"))
                    .spawn_scoped(scope, task)
                {
                    Ok(handle) => pending.push((index, handle)),
                    // Resource admission failure does not alter scientific
                    // failure semantics: this chunk uses the same serial path.
                    Err(_) => {
                        completed.push((index, chunk.iter().map(|p| curve(p.counts)).collect()))
                    }
                }
            }
            completed.push((0, first.iter().map(|p| curve(p.counts)).collect()));
            for (index, handle) in pending {
                match handle.join() {
                    Ok(values) => completed.push((index, values)),
                    Err(payload) => std::panic::resume_unwind(payload),
                }
            }
            completed.sort_unstable_by_key(|(index, _)| *index);
            completed
                .into_iter()
                .flat_map(|(_, values)| values)
                .collect()
        })
    };
    let first_likelihood = curve(base.first_category);
    let mut surface = vec![vec![0.0; 50]; 50];
    let mut best_index = [0, 0];
    let mut best = f64::NEG_INFINITY;
    for (i, row) in surface.iter_mut().enumerate() {
        let scale = 0.01 + (5.0 - 0.01) * i as f64 / 49.0;
        for (j, value) in row.iter_mut().enumerate() {
            let alpha = 10_f64.powf(-3.0 + (0.99_f64.log10() + 3.0) * j as f64 / 49.0);
            let p = response_prior(alpha, scale, step)?;
            let p0 = dot_dynamic(&first_likelihood, &p);
            *value = data
                .pairs
                .iter()
                .zip(&pair_likelihoods)
                .map(|(pair, l)| {
                    let q = dot_dynamic(l, &p) / (1.0 - p0);
                    if q > 0.0 {
                        pair.multiplicity as f64 / data.observations as f64 * q.ln()
                    } else {
                        0.0
                    }
                })
                .sum();
            if *value > best {
                best = *value;
                best_index = [i, j];
            }
        }
    }
    let alpha = 10_f64.powf(-3.0 + (0.99_f64.log10() + 3.0) * best_index[1] as f64 / 49.0);
    let scale = 0.01 + (5.0 - 0.01) * best_index[0] as f64 / 49.0;
    let reporting_prior = response_prior(alpha, scale, 0.1)?;
    let summaries = pair_likelihoods
        .iter()
        .map(|l| summarize(&shifts, &reporting_prior, l))
        .collect::<Result<Vec<_>, _>>()?;
    if surface.iter().flatten().any(|v| !v.is_finite()) {
        return Err("nonfinite legacy response surface".into());
    }
    Ok(Response {
        parameters1,
        parameters2,
        shifts,
        pair_likelihoods,
        first_likelihood,
        surface,
        best_index,
        alpha,
        scale,
        reporting_prior,
        summaries,
    })
}

/// Minimal legacy file workflow, using accepted strict ingestion and atomic
/// no-replace publication. Unsupported/nonconverged results publish nothing.
pub fn run_files(
    first: &std::path::Path,
    second: &std::path::Path,
    output: &std::path::Path,
    models: Option<[&std::path::Path; 2]>,
) -> Result<(), String> {
    use std::fs;
    let mut profile = crate::profile::Profile::new("legacy-files");
    let a = crate::read_sample(first).map_err(|e| e.to_string())?;
    profile.mark("sample1_read_hash_parse");
    let b = crate::read_sample(second).map_err(|e| e.to_string())?;
    profile.mark("sample2_read_hash_parse");
    let rows = crate::merge(&a, &b).map_err(|e| e.to_string())?;
    profile.mark("alignment");
    let data = Data::new(
        &rows
            .iter()
            .map(|r| [r.count1, r.count2])
            .collect::<Vec<_>>(),
    )?;
    let mut files: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    files.insert(
        "merged.tsv".into(),
        crate::encode(&rows).map_err(|e| e.to_string())?,
    );
    profile.mark("grouping_and_merged_table_encoding");
    let mut model_sources = Vec::new();
    let command = if let Some(paths) = models {
        let mut parameters = [[0.0; 2]; 2];
        for (i, path) in paths.iter().enumerate() {
            let bytes = fs::read(path).map_err(|e| e.to_string())?;
            let v: serde_json::Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
            if v["mode"] != "legacy-compatible"
                || v["upstream_revision"] != REVISION
                || v["fit"]["success"] != true
            {
                return Err(
                    "response requires successful legacy-compatible model.json files".into(),
                );
            }
            parameters[i] = serde_json::from_value(v["fit"]["evaluation"]["parameters"].clone())
                .map_err(|e| e.to_string())?;
            model_sources.push(serde_json::json!({"path":path.canonicalize().map_err(|e|e.to_string())?,"sha256":crate::sha256(&bytes)}));
        }
        let result = response(&data, parameters[0], parameters[1])?;
        profile.mark("response_compute");
        let mut table = String::from(
            "ntCDR3\tcount_1\tcount_2\ttail_P_s_le_zero\tmean\tmode\tlower\tmedian\tupper\n",
        );
        let mut selected = Vec::new();
        for row in &rows {
            let index = data
                .pairs
                .iter()
                .position(|p| p.counts == [row.count1, row.count2])
                .unwrap();
            let s = &result.summaries[index];
            let value = |v: Option<f64>| v.map_or("unavailable".into(), |v| format!("{v:.17}"));
            let line = format!(
                "{}\t{}\t{}\t{:.17}\t{:.17}\t{:.17}\t{}\t{}\t{}\n",
                row.nt,
                row.count1,
                row.count2,
                s.tail,
                s.mean,
                s.mode,
                value(s.lower),
                value(s.median),
                value(s.upper)
            );
            table.push_str(&line);
            if s.tail < 0.05 && s.median.is_some_and(|m| m > 0.0) {
                selected.push((s.tail, row.nt.clone(), line));
            }
        }
        selected.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut expanded = table.lines().next().unwrap().to_owned();
        expanded.push('\n');
        for (_, _, line) in selected {
            expanded.push_str(&line);
        }
        files.insert("response.tsv".into(), table.into_bytes());
        files.insert("expanded.tsv".into(), expanded.into_bytes());
        files.insert("response.json".into(), json_bytes(&result)?);
        "legacy-response"
    } else {
        let result = fit(&data, [-1.5, -4.0], 200)?;
        profile.mark("fit_compute");
        if !result.success {
            return Err(format!(
                "legacy fit nonconvergence: {}; C2={}",
                result.status, result.evaluation.c2
            ));
        }
        files.insert("diversity.json".into(),json_bytes(&serde_json::json!({"definition":"int(Nobs/(1-Pfirstfirst))","richness":result.evaluation.diversity}))?);
        let e = &result.evaluation;
        files.insert("diagnostics.json".into(),json_bytes(&serde_json::json!({"mean_nll":e.mean_nll,"C1":e.c1,"legacy_C2":e.c2,"first_category":e.first_category,"first_probability":e.first_probability,"explicit_zero_probability":e.explicit_zero_probability,"success":result.success,"status":result.status,"iterations":result.iterations}))?);
        files.insert(
            "parameters.tsv".into(),
            format!(
                "\tlabel\tvalue\n0\texponent\t{:.17}\n1\tlog10_minimum_frequency\t{:.17}\n",
                e.parameters[0], e.parameters[1]
            )
            .into_bytes(),
        );
        files.insert("model.json".into(),json_bytes(&serde_json::json!({"mode":"legacy-compatible","upstream_revision":REVISION,"fit":result}))?);
        "legacy-poisson"
    };
    profile.mark("scientific_output_encoding");
    let hashes: BTreeMap<_, _> = files
        .iter()
        .map(|(n, b)| (n.clone(), crate::sha256(b)))
        .collect();
    files.insert("provenance.json".into(),json_bytes(&serde_json::json!({
        "mode":"legacy-compatible","command":command,"version":crate::VERSION,"upstream_revision":REVISION,
        "inputs":[{"path":a.path,"sha256":a.sha256,"depth":a.depth},{"path":b.path,"sha256":b.sha256,"depth":b.depth}],
        "models":model_sources,"observed_union":rows.len(),"data":data,"validation_policy":crate::POLICY,"outputs_sha256":hashes,
        "settings":{"grid_points":POINTS,"log_cutoff_base":10,"prior_normalization":"legacy trapezoid","conditioning":"first observed axis categories","constraint":"legacy defective C2 equality; C1 diagnostic","optimizer":"native equality SQP/BFGS, central derivative 1e-5, 200 iterations, equality <1e-6, stationarity <1e-7","initial":[-1.5,-4.0],"parameter_bounds":"none","response":"501 aligned shifts, 50x50 hyperparameter grid; nominal .1 reporting prior","expansion_filter":{"tail_strictly_below":0.05,"median_strictly_above":0.0}},
        "protocol_sha256":"7fc253760a4a1a24c6308685cf9f60461526d60507c241c900318eddedc09b11",
        "implementation_settings_sha256":"8679951fed3784a10180ae9664cafe3b64bee5628b1e2ce299510025e0ddb868",
        "limitations":["known first-category and C2 defects reproduced","optimizer trajectory differs from SLSQP","no calibrated frequentist p-values","no legacy contraction selection","no NB/NBP fit or generator","local qualification only"]
    }))?);
    profile.mark("output_hashes_and_provenance");
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(std::path::Path::new("."))
        .canonicalize()
        .map_err(|e| e.to_string())?;
    let name = output
        .file_name()
        .ok_or("output directory needs a new final name")?;
    let destination = parent.join(name);
    if fs::symlink_metadata(&destination).is_ok() {
        return Err("refusing existing legacy output directory".into());
    }
    let stage = tempfile::Builder::new()
        .prefix(".clonoshift-legacy-")
        .tempdir_in(parent)
        .map_err(|e| e.to_string())?;
    for (name, bytes) in files {
        crate::write_synced(&stage.path().join(name), &bytes).map_err(|e| e.to_string())?;
    }
    fs::File::open(stage.path())
        .and_then(|f| f.sync_all())
        .map_err(|e| e.to_string())?;
    crate::publish(stage.path(), &destination).map_err(|e| e.to_string())?;
    profile.mark("synchronized_atomic_publication");
    Ok(())
}
fn json_bytes(v: &impl Serialize) -> Result<Vec<u8>, String> {
    let mut bytes = serde_json::to_vec_pretty(v).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    Ok(bytes)
}
