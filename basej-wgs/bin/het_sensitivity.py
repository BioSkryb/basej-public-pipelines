#!/usr/bin/env python3
"""Theoretical het-SNP sensitivity, ported from Picard's TheoreticalSensitivity.

basej-wgs gets HET_SNP_SENSITIVITY / HET_SNP_Q for free: Sentieon's WgsMetricsAlgo writes
them into its metrics file and the pipeline just reads the columns. bskryb-qc has no such
model, so the two fields were the only values in the shared 97-field wgsqc_summary schema
that the Rust producer left null. This reimplements the model over the two histograms
bskryb-qc now emits (quality_yield.QUALITY_HISTOGRAM, coverage.DEPTH_HISTOGRAM).

Faithful to picard/analysis/TheoreticalSensitivity.java as called from WgsMetrics.java:

    depthDoubleArray = normalizeHistogram(depthHistogram)     # no trimDistribution
    baseQDoubleArray = normalizeHistogram(baseQHistogram)
    HET_SNP_SENSITIVITY = hetSNPSensitivity(depth, baseQ, SAMPLE_SIZE, LOG_ODDS_THRESHOLD=3.0)
    HET_SNP_Q           = QualityUtil.getPhredScoreFromErrorProbability(1 - sensitivity)

Not bit-exact with Picard, and not expected to be: the estimator is Monte Carlo over a
stochastic-acceptance roulette wheel seeded with java.util.Random(51), and Java's RNG
stream cannot be reproduced from numpy. Validated instead against Sentieon's published
values -- feeding Sentieon's own depth histogram for ResolveOMEv2.0-384-WGS-B4 reproduces
0.96483 against its reported 0.96628 (-0.15%), with HET_SNP_Q identical at 15. Monte Carlo
noise is not the limiting factor there: across 8 seeds the spread is 1e-4.

One deliberate deviation: Picard feeds its *unfiltered* depth histogram, whereas
coverage.DEPTH_HISTOGRAM is the high-quality (post MAPQ/baseQ/overlap) depth that also
produces our MEAN_COVERAGE and PCT_*X. Using it keeps every coverage-derived field on one
consistent depth definition, and biases sensitivity slightly low (conservative) rather than
high -- scaling depth up by B4's PCT_EXC_TOTAL overshoots Sentieon's value.
"""
import math

import numpy as np

# Picard TheoreticalSensitivity constants.
MAX_CONSIDERED_DEPTH = 1000  # MAX_CONSIDERED_DEPTH_HET_SENS
RANDOM_SEED = 51
LOG_ODDS_THRESHOLD = 3.0  # WgsMetrics.LOG_ODDS_THRESHOLD
SAMPLE_SIZE = 10000  # CollectWgsMetrics THEORETICAL_SENSITIVITY_SAMPLE_SIZE default


def _normalized(hist):
    """Picard normalizeHistogram: divide by the sum. None when there is nothing to divide."""
    if hist is None:
        return None
    arr = np.asarray(hist, dtype=np.float64)
    if arr.size == 0 or not np.isfinite(arr).all() or arr.sum() <= 0:
        return None
    return arr / arr.sum()


def het_snp_sensitivity(depth_histogram, quality_histogram,
                        sample_size=SAMPLE_SIZE,
                        log_odds_threshold=LOG_ODDS_THRESHOLD,
                        seed=RANDOM_SEED):
    """Probability of calling a het SNP, integrated over the observed depth distribution.

    Returns None if either histogram is missing or empty, so that "not measured" stays
    distinguishable from "measured and low" -- the same contract qc_scoring.qc_numeric keeps.
    """
    depth = _normalized(depth_histogram)
    quality = _normalized(quality_histogram)
    if depth is None or quality is None:
        return None

    n_depths = min(len(depth), MAX_CONSIDERED_DEPTH + 1)
    depth = depth[:n_depths]

    # qualitySums[:, m] is a sample of sums of m quality draws, so column 0 is identically
    # zero (a het site with zero alt reads contributes no quality evidence).
    rng = np.random.default_rng(seed)
    draws = rng.choice(len(quality), size=(sample_size, n_depths), p=quality)
    quality_sums = np.zeros((sample_size, n_depths), dtype=np.int64)
    if n_depths > 1:
        np.cumsum(draws[:, : n_depths - 1], axis=1, out=quality_sums[:, 1:])

    # A SNP is called when the summed alt-base qualities clear this threshold. Picard's
    # `LOG_10` local is log10(2), not log10(10) -- the 2 is the het allele fraction.
    thresholds = 10.0 * (np.arange(n_depths) * math.log10(2.0) + log_odds_threshold)

    # exceed[m, n] = fraction of the m-summand sample at or above thresholds[n]. Picard's
    # sorted scan advances while threshold > value, making the boundary inclusive;
    # searchsorted(side='left') counts strictly-smaller values, which is the same split.
    exceed = np.empty((n_depths, n_depths), dtype=np.float64)
    for m in range(n_depths):
        ordered = np.sort(quality_sums[:, m])
        below = np.searchsorted(ordered, thresholds, side="left")
        exceed[m] = (sample_size - below) / sample_size

    # alt_depth[n, m] = C(n, m) * 0.5**n -- the chance a het site with depth n yields m alt
    # reads. Built with Picard's nCm = (n-1)C(m-1) * (n/m) recurrence to avoid overflowing
    # binomial coefficients. Entries with m > n stay exactly zero, which is what makes the
    # m <= n restriction of the double sum implicit in the einsum below.
    alt_depth = np.zeros((n_depths, n_depths), dtype=np.float64)
    alt_depth[0, 0] = 1.0
    for n in range(1, n_depths):
        alt_depth[n, 0] = 0.5 ** n
        if n > 1:
            m = np.arange(1, n)
            alt_depth[n, 1:n] = (n * 0.5 / m) * alt_depth[n - 1, 0 : n - 1]
        alt_depth[n, n] = alt_depth[n, 0]

    # sum over n of depth[n] * sum over m <= n of alt_depth[n, m] * exceed[m, n]
    return float(depth @ np.einsum("nm,mn->n", alt_depth, exceed))


def het_snp_q(sensitivity):
    """Picard QualityUtil.getPhredScoreFromErrorProbability, which rounds to an int."""
    if sensitivity is None:
        return None
    error_rate = 1.0 - sensitivity
    if error_rate <= 0.0:
        return None  # saturated model; a phred of +inf is not a value worth publishing
    return int(round(-10.0 * math.log10(error_rate)))
