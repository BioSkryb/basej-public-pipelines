#!/usr/bin/env python3
"""Apply a chimera-class ruleset (chimera_classes.yaml) to a bskryb-qc rustqc JSON's
`chimera_cube`, producing mutually-exclusive per-class read counts and fractions.

This is the modular classification layer: the Rust binary emits the cube (feature
substrate); named mechanism classes live here as declarative rules, so new/experimental
classes are a config edit, not a re-run. Reused by the pipeline parquet writers.

    from chimera_classify import classify_cube, load_ruleset
    classes = classify_cube(rustqc_json["chimera_cube"], load_ruleset("chimera_classes.yaml"))
    # -> {"proper": {"count":..,"pct":..,"chimeric":False}, "interchromosomal": {...}, ...}
"""
import json
import os
import sys

_DEFAULT_RULES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "chimera_classes.yaml")

# Embedded v1 ruleset — mirrors chimera_classes.yaml so the pipeline parquet writers can
# classify without a pyyaml dependency in the metrics container. The YAML remains the
# canonical, editable source (used by this tool for exploration); keep them in sync.
_V1_RULESET = {
    "version": "v1",
    "classes": [
        {"name": "proper", "chimeric": False,
         "match": {"locus": ["same_contig"], "orient": ["proper_FR"],
                   "dist": ["lt_1kb", "1_10kb", "10_100kb"], "sa": ["no_SA"]}},
        {"name": "interchromosomal", "match": {"locus": ["diff_contig"]}},
        {"name": "local_inverted_hairpin",
         "match": {"locus": ["same_contig"], "orient": ["FF", "RR"], "dist": ["lt_1kb"]}},
        {"name": "inverted_nonlocal",
         "match": {"locus": ["same_contig"], "orient": ["FF", "RR"],
                   "dist": ["1_10kb", "10_100kb", "100kb_1Mb", "gt_1Mb"]}},
        {"name": "everted_tandem",
         "match": {"locus": ["same_contig"], "orient": ["everted"]}},
        {"name": "large_insert",
         "match": {"locus": ["same_contig"], "orient": ["proper_FR"], "dist": ["100kb_1Mb", "gt_1Mb"]}},
        {"name": "split_read_only",
         "match": {"locus": ["same_contig"], "orient": ["proper_FR"],
                   "dist": ["lt_1kb", "1_10kb", "10_100kb"], "sa": ["has_SA"]}},
        {"name": "other", "match": {}},
    ],
}


def load_ruleset(path=_DEFAULT_RULES):
    """Load the ruleset from YAML if pyyaml + the file are available; otherwise fall back
    to the embedded v1 ruleset (keeps the pipeline independent of pyyaml)."""
    try:
        import yaml
        if path and os.path.exists(path):
            with open(path) as fh:
                return yaml.safe_load(fh)
    except ImportError:
        pass
    return _V1_RULESET


def _match(cell, m):
    """A cell matches a rule when, for every axis the rule names, the cell's value is listed."""
    for axis, allowed in m.items():
        if cell.get(axis) not in allowed:
            return False
    return True


def classify_cube(cube, ruleset):
    """Assign each cube cell to the first matching class (priority order). Returns a dict
    of class_name -> {count, pct, chimeric}. pct is fraction of reads_aligned_in_pairs."""
    n_pairs = cube.get("reads_aligned_in_pairs", 0) or 0
    rules = ruleset["classes"]
    out = {r["name"]: {"count": 0, "chimeric": bool(r.get("chimeric", True))} for r in rules}
    for cell in cube.get("cells", []):
        cnt = cell.get("count", 0)
        for r in rules:
            if _match(cell, r.get("match", {})):
                out[r["name"]]["count"] += cnt
                break
    for name, rec in out.items():
        rec["pct"] = (rec["count"] / n_pairs) if n_pairs else 0.0
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: chimera_classify.py <rustqc.json> [ruleset.yaml]")
    with open(sys.argv[1]) as fh:
        j = json.load(fh)
    rules = load_ruleset(sys.argv[2] if len(sys.argv) > 2 else _DEFAULT_RULES)
    cube = j["chimera_cube"]
    res = classify_cube(cube, rules)
    n = cube.get("reads_aligned_in_pairs", 0)
    print(f"ruleset {rules.get('version')}   reads_aligned_in_pairs={n:,}")
    chim_total = 0
    for name, rec in res.items():
        tag = "chimeric" if rec["chimeric"] else "proper  "
        if rec["chimeric"]:
            chim_total += rec["count"]
        print(f"  {tag}  {name:<24} {rec['count']:>10,}  {rec['pct']*100:6.3f}%")
    print(f"  ---- total chimeric = {chim_total:,}  ({chim_total/max(n,1)*100:.3f}% = PCT_CHIMERAS)")


if __name__ == "__main__":
    main()
