#!/usr/bin/env python3
"""
Run MuSiCaL COSMIC refitting on a SBS96 mutation count matrix.

Loads the COSMIC v3.2 SBS WGS reference signatures directly from MuSiCaL's
built-in catalog (COSMIC_v3p2_SBS_WGS) and refits each sample via
musical.refit.refit() using the likelihood_bidirectional sparse-NNLS method.

WES and WGS both use the same catalog; no WES-specific catalog is available
in MuSiCaL (same limitation as the SigDyn pipeline).

Usage:
    run_musical_analysis.py <sbs96_matrix> <cohort_id> <genome_build> <exome> <output_prefix>

Arguments:
    sbs96_matrix   Tab-separated SBS96 matrix (96 rows x N samples) from
                   SigProfilerMatrixGenerator (MutationType column + one col per sample)
    cohort_id      Cohort identifier used in output filenames
    genome_build   GRCh37 or GRCh38 (informational only)
    exome          true/false (informational; same catalog used for both modes)
    output_prefix  Prefix for all output files
"""
import sys
import pandas as pd
import numpy as np


COSMIC_CATALOG_NAME = "COSMIC_v3p2_SBS_WGS"


def _cos_sim(a, b):
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denom) if denom > 0.0 else 0.0


def main():
    if len(sys.argv) < 6:
        print(
            "Usage: run_musical_analysis.py "
            "<sbs96_matrix> <cohort_id> <genome_build> <exome> <output_prefix>",
            file=sys.stderr,
        )
        sys.exit(1)

    matrix_file  = sys.argv[1]
    cohort_id    = sys.argv[2]
    genome_build = sys.argv[3]
    exome        = sys.argv[4].lower() in ("true", "1", "yes")
    out_prefix   = sys.argv[5]

    import musical
    import musical.refit as refit_mod
    from musical.catalog import load_catalog

    # ── Load SBS96 count matrix ───────────────────────────────────────────────
    print(f"Loading SBS96 matrix: {matrix_file}", flush=True)
    X = pd.read_csv(matrix_file, sep="\t", index_col=0)
    print(f"Matrix shape: {X.shape[0]} mutation types × {X.shape[1]} samples", flush=True)

    if X.shape[1] == 0:
        print("ERROR: no samples found in matrix", file=sys.stderr)
        sys.exit(1)

    # ── Load COSMIC reference catalog from MuSiCaL built-in ──────────────────
    # The pre-built TSV exported from MutationalPatterns has numeric row indices
    # instead of mutation type names (get_known_signatures() returns a matrix
    # with NULL rownames), so we load directly from MuSiCaL where features are
    # properly labelled (e.g. "A[C>A]A").
    # No WES-specific catalog is available in MuSiCaL; WES mode falls back to WGS.
    mode_label = "WES" if exome else "WGS"
    print(f"Loading COSMIC v3.2 catalog ({mode_label} mode) from MuSiCaL: {COSMIC_CATALOG_NAME}", flush=True)
    cat = load_catalog(COSMIC_CATALOG_NAME)
    W = pd.DataFrame(cat.W, index=cat.features, columns=cat.signatures)
    print(f"Catalog: {W.shape[0]} mutation types × {W.shape[1]} signatures", flush=True)

    # ── Align mutation type order ─────────────────────────────────────────────
    common = X.index.intersection(W.index)
    if len(common) < 90:
        print(
            f"WARNING: only {len(common)}/96 mutation types matched — "
            "verify input is a SigProfilerMatrixGenerator SBS96 matrix",
            file=sys.stderr,
        )
    X_aln = X.loc[common].astype(float)
    W_aln = W.loc[common].astype(float)

    # ── MuSiCaL likelihood-bidirectional refitting ────────────────────────────
    print(
        f"Running MuSiCaL refitting  "
        f"({X_aln.shape[1]} samples, {W_aln.shape[1]} COSMIC sigs, {genome_build})...",
        flush=True,
    )
    H_raw, _model = refit_mod.refit(
        X_aln,
        W_aln,
        method="likelihood_bidirectional",
    )
    # H_raw: K signatures × N samples (DataFrame or ndarray)
    H_np = H_raw.values if hasattr(H_raw, "values") else np.array(H_raw)

    # ── Export signature activities (samples × signatures) ────────────────────
    sig_names   = W_aln.columns.tolist()
    sample_names = X_aln.columns.tolist()
    activities  = pd.DataFrame(H_np.T, index=sample_names, columns=sig_names)
    activities.index.name = "Samples"
    act_file = f"{out_prefix}_musical_activities.tsv"
    activities.to_csv(act_file, sep="\t")
    print(f"Saved activities  → {act_file}", flush=True)

    # ── Per-sample cosine similarity ──────────────────────────────────────────
    recon    = W_aln.values @ H_np        # 96 × N
    orig     = X_aln.values
    cos_sims = np.array([_cos_sim(orig[:, i], recon[:, i]) for i in range(orig.shape[1])])
    cos_df   = pd.DataFrame({"Sample": sample_names, "CosineSimilarity": cos_sims})
    cos_file = f"{out_prefix}_musical_cosine_similarities.tsv"
    cos_df.to_csv(cos_file, sep="\t", index=False)
    print(f"Saved cosine sims → {cos_file}", flush=True)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_active = int((activities > 0).any(axis=0).sum())
    print(
        f"Active COSMIC signatures: {n_active}/{activities.shape[1]}   "
        f"Mean cosine similarity: {cos_sims.mean():.4f}",
        flush=True,
    )
    print("MuSiCaL analysis complete.", flush=True)


if __name__ == "__main__":
    main()
