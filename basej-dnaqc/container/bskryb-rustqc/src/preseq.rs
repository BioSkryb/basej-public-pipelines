// preseq gc_extrap reimplementation (library/genome-coverage complexity).
//
// Faithful Rust port of the deterministic core of preseq gc_extrap:
//   - continued_fraction.cpp  (QD algorithm, Euler evaluation, stability search)
//   - common.cpp              (Lanczos gamma, Heck interpolation, extrapolate_curve,
//                              extrap_single_estimate, bootstrap + median)
// Original: Daley & Smith, preseq (smithlabcode/preseq, GPLv3+).
//
// Difference from upstream: bit-exact reproduction is impossible (upstream depends on
// GSL/libstdc++ RNG streams for both the probabilistic binning and the bootstrap). We
// keep the SAME methodology (100-resample bootstrap median) with a fixed-seed Rust RNG,
// so the result is deterministic within this tool and close to upstream. AI-assisted;
// validated against the preseq oracle within tolerance.

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rand_distr::{Binomial, Distribution};

const MIN_ALLOWED_DEGREE: usize = 4;
const SEARCH_MAX_VAL: f64 = 100.0;
const SEARCH_STEP_SIZE: f64 = 0.05;

/// Σ_j j·hist[j]  (total coverage events / observed sample size in bin units).
fn counts_from_hist(hist: &[f64]) -> f64 {
    hist.iter().enumerate().map(|(i, &c)| i as f64 * c).sum()
}

pub fn good_toulmin_2x(hist: &[f64]) -> f64 {
    let mut s = 0.0;
    for (i, &c) in hist.iter().enumerate() {
        s += (-1.0f64).powi((i + 1) as i32) * c;
    }
    s
}

/// Lanczos log-gamma; matches preseq common.cpp `factorial` (log((x-1)!)).
fn factorial(mut x: f64) -> f64 {
    const LOG_ROOT_TWO_PI: f64 = 0.9189385332046727;
    const EULER: f64 = 2.71828182845904523536028747135;
    const LANCZOS: [f64; 9] = [
        0.99999999999980993227684700473478,
        676.520368121885098567009190444019,
        -1259.13921672240287047156078755283,
        771.3234287776530788486528258894,
        -176.61502916214059906584551354,
        12.507343278686904814458936853,
        -0.13857109526572011689554707,
        9.984369578019570859563e-6,
        1.50563273514931155834e-7,
    ];
    x -= 1.0;
    let mut ag = LANCZOS[0];
    for (k, &l) in LANCZOS.iter().enumerate().skip(1) {
        ag += l / (x + k as f64);
    }
    let term1 = (x + 0.5) * ((x + 7.5) / EULER).ln();
    let term2 = LOG_ROOT_TWO_PI + ag.ln();
    term1 + (term2 - 7.0)
}

/// Heck 1975 exact rarefaction: expected distinct in subsample n of N with S distinct.
fn interpolate_distinct(hist: &[f64], n_total: usize, s: usize, sub: usize) -> f64 {
    let n = n_total as f64;
    let nn = sub as f64;
    let denom = factorial(n + 1.0) - factorial(nn + 1.0) - factorial(n - nn + 1.0);
    let mut numer_sum = 0.0;
    for (i, &h) in hist.iter().enumerate() {
        if i == 0 {
            continue;
        }
        if n_total < i + sub {
            continue;
        }
        let x = factorial(n - i as f64 + 1.0) - factorial(nn + 1.0)
            - factorial(n - i as f64 - nn + 1.0);
        numer_sum += (x - denom).exp() * h;
    }
    s as f64 - numer_sum
}

// ------------------------- Continued fraction -------------------------

#[derive(Clone, Default)]
struct ContinuedFraction {
    ps_coeffs: Vec<f64>,
    cf_coeffs: Vec<f64>,
    offset_coeffs: Vec<f64>,
    diagonal_idx: i32,
    degree: usize,
}

fn quotdiff_algorithm(ps: &[f64]) -> Vec<f64> {
    let depth = ps.len();
    if depth == 0 {
        return vec![];
    }
    let mut q = vec![vec![0.0f64; depth + 1]; depth];
    for j in 0..depth.saturating_sub(1) {
        q[1][j] = ps[j + 1] / ps[j];
    }
    let mut e = vec![vec![0.0f64; depth + 1]; depth];
    for j in 0..depth.saturating_sub(1) {
        e[1][j] = q[1][j + 1] - q[1][j] + e[0][j + 1];
    }
    for i in 2..depth {
        for j in 0..depth {
            q[i][j] = q[i - 1][j + 1] * e[i - 1][j + 1] / e[i - 1][j];
        }
        for j in 0..depth {
            e[i][j] = q[i][j + 1] - q[i][j] + e[i - 1][j + 1];
        }
    }
    let mut cf = vec![0.0f64; depth];
    cf[0] = ps[0];
    for i in 1..depth {
        cf[i] = if i % 2 == 0 { -e[i / 2][0] } else { -q[(i + 1) / 2][0] };
    }
    cf
}

fn get_rescale_value(num: f64, den: f64) -> f64 {
    const TOL: f64 = 1e-20;
    let rescale = num.abs() + den.abs();
    if rescale > 1.0 / TOL {
        1.0 / rescale
    } else if rescale < TOL {
        1.0 / rescale
    } else {
        1.0
    }
}

fn evaluate_on_diagonal(cf_coeffs: &[f64], val: f64, depth: usize) -> f64 {
    let mut current_num;
    let mut prev_num1 = cf_coeffs[0];
    let mut prev_num2 = 0.0;
    let mut current_denom;
    let mut prev_denom1 = 1.0;
    let mut prev_denom2 = 1.0;
    let lim = cf_coeffs.len().min(depth);
    // if lim <= 1, no iterations; return cf_coeffs[0]/1
    let mut cn = prev_num1;
    let mut cd = prev_denom1;
    for i in 1..lim {
        current_num = prev_num1 + cf_coeffs[i] * val * prev_num2;
        current_denom = prev_denom1 + cf_coeffs[i] * val * prev_denom2;
        prev_num2 = prev_num1;
        prev_num1 = current_num;
        prev_denom2 = prev_denom1;
        prev_denom1 = current_denom;
        let rescale = get_rescale_value(current_num, current_denom);
        prev_num1 *= rescale;
        prev_num2 *= rescale;
        prev_denom1 *= rescale;
        prev_denom2 *= rescale;
        cn = prev_num1;
        cd = prev_denom1;
    }
    cn / cd
}

fn evaluate_power_series(ps: &[f64], val: f64) -> f64 {
    ps.iter().enumerate().map(|(i, &c)| c * val.powi(i as i32)).sum()
}

impl ContinuedFraction {
    fn new(ps: Vec<f64>, di: i32, dg: usize) -> Self {
        let mut cf = ContinuedFraction {
            ps_coeffs: ps,
            cf_coeffs: vec![],
            offset_coeffs: vec![],
            diagonal_idx: di,
            degree: dg,
        };
        if di == 0 {
            cf.cf_coeffs = quotdiff_algorithm(&cf.ps_coeffs);
        } else if di > 0 {
            let offset = di as usize;
            let high: Vec<f64> = cf.ps_coeffs[offset..].to_vec();
            cf.cf_coeffs = quotdiff_algorithm(&high);
            cf.offset_coeffs = cf.ps_coeffs[..offset].to_vec();
        } else {
            let offset = (-di) as usize;
            let n = cf.ps_coeffs.len();
            let mut recip = vec![0.0f64; n];
            recip[0] = 1.0 / cf.ps_coeffs[0];
            for i in 1..n {
                let mut x = 0.0;
                for j in 0..i {
                    x += cf.ps_coeffs[i - j] * recip[j];
                }
                recip[i] = -x / cf.ps_coeffs[0];
            }
            let high: Vec<f64> = recip[offset..].to_vec();
            cf.cf_coeffs = quotdiff_algorithm(&high);
            cf.offset_coeffs = recip[..offset].to_vec();
        }
        cf
    }

    fn is_valid(&self) -> bool {
        !self.cf_coeffs.is_empty()
    }

    fn eval(&self, val: f64) -> f64 {
        if self.diagonal_idx > 0 {
            let cf_part = evaluate_on_diagonal(&self.cf_coeffs, val, self.degree - self.offset_coeffs.len());
            let ps_part = evaluate_power_series(&self.offset_coeffs, val);
            ps_part + val.powi(self.offset_coeffs.len() as i32) * cf_part
        } else if self.diagonal_idx < 0 {
            let cf_part = evaluate_on_diagonal(&self.cf_coeffs, val, self.degree - self.offset_coeffs.len());
            let ps_part = evaluate_power_series(&self.offset_coeffs, val);
            1.0 / (ps_part + val.powi(self.offset_coeffs.len() as i32) * cf_part)
        } else {
            evaluate_on_diagonal(&self.cf_coeffs, val, self.degree)
        }
    }

    fn truncate(&self, n_terms: usize) -> Option<ContinuedFraction> {
        if self.degree < n_terms {
            None
        } else {
            let mut c = self.clone();
            c.ps_coeffs.truncate(n_terms);
            let new_cf_len = n_terms - c.offset_coeffs.len();
            c.cf_coeffs.truncate(new_cf_len);
            c.degree = n_terms;
            Some(c)
        }
    }

    fn extrapolate_distinct(&self, max_value: f64, step: f64) -> Vec<f64> {
        let mut est = vec![0.0];
        let mut t = step;
        while t <= max_value {
            est.push(t * self.eval(t));
            t += step;
        }
        est
    }
}

fn check_stability(est: &[f64]) -> bool {
    if est.is_empty() {
        return false;
    }
    for &e in est {
        if !e.is_finite() || e < 0.0 {
            return false;
        }
    }
    for i in 1..est.len() {
        if est[i] < est[i - 1] {
            return false;
        }
    }
    for i in 2..est.len() {
        if est[i - 1] - est[i - 2] < est[i] - est[i - 1] {
            return false;
        }
    }
    true
}

/// optimal_cont_frac_distinct: find largest stable CF degree.
fn optimal_cont_frac_distinct(hist: &[f64], diagonal: i32, max_terms: usize) -> ContinuedFraction {
    if max_terms >= hist.len() {
        return ContinuedFraction::default();
    }
    let mut ps: Vec<f64> = Vec::new();
    for j in 1..=max_terms {
        ps.push(hist[j] * (-1.0f64).powi((j + 1) as i32));
    }
    let full = ContinuedFraction::new(ps, diagonal, max_terms);
    if (3..=6).contains(&max_terms) {
        let est = full.extrapolate_distinct(SEARCH_MAX_VAL, SEARCH_STEP_SIZE);
        if check_stability(&est) {
            return full;
        }
    } else {
        let start = 7 + (max_terms % 2 == 0) as usize;
        let mut i = start;
        while i <= max_terms {
            if let Some(trunc) = full.truncate(i) {
                let est = trunc.extrapolate_distinct(SEARCH_MAX_VAL, SEARCH_STEP_SIZE);
                if check_stability(&est) {
                    return trunc;
                }
            }
            i += 2;
        }
    }
    ContinuedFraction::default()
}

/// extrapolate_curve: estimates = initial_distinct + fold·CF(fold).
fn extrapolate_curve(
    cf: &ContinuedFraction,
    initial_distinct: f64,
    vals_sum: f64,
    initial_sample: f64,
    step: f64,
    max_sample: f64,
    est: &mut Vec<f64>,
) {
    let mut curr = initial_sample;
    while curr < max_sample {
        let fold = (curr - vals_sum) / vals_sum;
        est.push(initial_distinct + fold * cf.eval(fold));
        curr += step;
    }
}

/// One deterministic estimate curve (the -Q path). Returns None if no stable CF.
fn single_estimate(
    hist: &[f64],
    mut max_terms: usize,
    diagonal: i32,
    step: f64,
    max_extrap: f64,
) -> Option<Vec<f64>> {
    let vals_sum = counts_from_hist(hist);
    let initial_distinct: f64 = hist.iter().sum();
    let mut yield_est: Vec<f64> = Vec::new();

    let upper = vals_sum as usize;
    let step_us = step as usize;
    let mut sample = step as usize;
    while sample < upper {
        yield_est.push(interpolate_distinct(hist, upper, initial_distinct as usize, sample));
        sample += step_us;
    }

    let mut first_zero = 1usize;
    while first_zero < hist.len() && hist[first_zero] > 0.0 {
        first_zero += 1;
    }
    max_terms = max_terms.min(first_zero - 1);
    if max_terms % 2 == 1 {
        max_terms -= 1;
    }
    if max_terms < MIN_ALLOWED_DEGREE {
        return None;
    }

    let cf = optimal_cont_frac_distinct(hist, diagonal, max_terms);
    if !cf.is_valid() {
        return None;
    }
    extrapolate_curve(&cf, initial_distinct, vals_sum, sample as f64, step, max_extrap, &mut yield_est);
    Some(yield_est)
}

/// Multinomial draw of `n` into buckets with weights `probs` (proportional).
fn multinomial(rng: &mut ChaCha8Rng, probs: &[f64], n: u64) -> Vec<u64> {
    let total: f64 = probs.iter().sum();
    let mut out = vec![0u64; probs.len()];
    let mut remaining = n;
    let mut remaining_p = total;
    for (i, &p) in probs.iter().enumerate() {
        if remaining == 0 || remaining_p <= 0.0 {
            break;
        }
        if i == probs.len() - 1 {
            out[i] = remaining;
            break;
        }
        let prob = (p / remaining_p).clamp(0.0, 1.0);
        let draw = Binomial::new(remaining, prob).unwrap().sample(rng);
        out[i] = draw;
        remaining -= draw;
        remaining_p -= p;
    }
    out
}

/// Bootstrap: resample histogram, run single_estimate per resample, keep stable ones.
fn bootstrap(
    hist: &[f64],
    n_bootstraps: usize,
    max_terms: usize,
    diagonal: i32,
    step: f64,
    max_extrap: f64,
    seed: u64,
    max_iter: usize,
) -> Vec<Vec<f64>> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    // distinct (nonzero) entries
    let mut idx: Vec<usize> = Vec::new();
    let mut vals: Vec<f64> = Vec::new();
    for (i, &c) in hist.iter().enumerate() {
        if c > 0.0 {
            idx.push(i);
            vals.push(c);
        }
    }
    let distinct: u64 = vals.iter().sum::<f64>() as u64;
    let mut out: Vec<Vec<f64>> = Vec::new();
    let mut iter = 0;
    while iter < max_iter && out.len() < n_bootstraps {
        iter += 1;
        let sample = multinomial(&mut rng, &vals, distinct);
        let mut boot = vec![0.0f64; idx[idx.len() - 1] + 1];
        for (k, &i) in idx.iter().enumerate() {
            boot[i] = sample[k] as f64;
        }
        while boot.len() > 1 && *boot.last().unwrap() == 0.0 {
            boot.pop();
        }
        if let Some(est) = single_estimate(&boot, max_terms, diagonal, step, max_extrap) {
            if check_stability(&est) {
                out.push(est);
            }
        }
    }
    out
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// gc_extrap estimate: returns the extrapolated COVERED BASES at the max step
/// (bootstrap median), i.e. the pipeline's `preseq_count`. None if it can't estimate.
pub fn gc_extrap_estimate(
    coverage_hist: &[f64],
    bin_size: usize,
    base_step_size: f64,
    max_extrap: f64,
    n_bootstraps: usize,
    seed: u64,
) -> Option<f64> {
    if coverage_hist.len() < 2 {
        return None;
    }
    // saturation guard
    if good_toulmin_2x(coverage_hist) < 0.0 {
        return None;
    }
    // ensure enough terms before first zero
    let mut first_zero = 1usize;
    while first_zero < coverage_hist.len() && coverage_hist[first_zero] > 0.0 {
        first_zero += 1;
    }
    let orig_max_terms = 100usize.min(first_zero - 1);
    if orig_max_terms < MIN_ALLOWED_DEGREE {
        return None;
    }
    let bin_step = base_step_size / bin_size as f64;
    let max_ex = max_extrap / bin_size as f64;

    let boot = bootstrap(
        coverage_hist,
        n_bootstraps,
        orig_max_terms,
        0,
        bin_step,
        max_ex,
        seed,
        10 * n_bootstraps,
    );
    if std::env::var("PRESEQ_DEBUG").is_ok() {
        let single = single_estimate(coverage_hist, orig_max_terms, 0, bin_step, max_ex);
        let single_tail = single.as_ref().and_then(|v| v.last().copied()).unwrap_or(-1.0);
        eprintln!(
            "[preseq debug] orig_max_terms={orig_max_terms} bin_step={bin_step} max_ex={max_ex} boots_ok={} single_tail_bins={single_tail} single_tail_bases={}",
            boot.len(),
            single_tail * bin_size as f64
        );
    }
    if boot.is_empty() {
        return None;
    }
    // median of the last column (the max extrapolation point), × bin_size
    let min_len = boot.iter().map(|b| b.len()).min().unwrap_or(0);
    if min_len == 0 {
        return None;
    }
    let last = min_len - 1;
    let col: Vec<f64> = boot.iter().map(|b| b[last]).collect();
    Some(median(col) * bin_size as f64)
}
