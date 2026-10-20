#!/usr/bin/env bash
# Build all x86 basej-dnaqc containers locally, one log per image + a summary.
# Image tags match the `container` directives in ../nextflow.config. fastp,
# bedtools and MultiQC are pulled from quay.io/biocontainers and are not built here.
# Failures are recorded as FAILED in _build_logs/SUMMARY.log (the script itself
# does not abort, so every image gets a build attempt and a log).
set -u
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="$BASE/_build_logs"
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/SUMMARY.log"
: > "$SUMMARY"
# folder:image_tag  (x86 only)
BUILDS=(
  "seqkit:basejumper_seqkit-2.13.0"
  "samtools:basejumper_samtools-1.23.1"
  "minibwa-rustqc:basejumper_minibwa-rustqc_0.1.5"
  "bskryb-rustqc:basejumper_bskryb-rustqc-x86_0.1.5"
  "preseq_bam2mr:basejumper_preseq_bam2mr_0.1"
  "ginkgo:basejumper_ginkgo_0.3.1"
  "ginko_parser:basejumper_ginko_parser_0.2.1"
  "custom_parabricks-metrics:basejumper_custom_parabricks-metrics_1.0.3"
  "custom_r_qcplots:basejumper_custom_r_qcplots_0.3.0"
)
echo "BUILD STARTED: $(date)" | tee -a "$SUMMARY"
for entry in "${BUILDS[@]}"; do
  folder="${entry%%:*}"
  tag="${entry##*:}"
  log="$LOGDIR/${folder}.log"
  echo "[$(date +%H:%M:%S)] BUILDING $folder -> $tag (log: $log)" | tee -a "$SUMMARY"
  if docker build --platform linux/amd64 -t "$tag" "$BASE/$folder" > "$log" 2>&1; then
    echo "[$(date +%H:%M:%S)] OK      $folder -> $tag" | tee -a "$SUMMARY"
  else
    echo "[$(date +%H:%M:%S)] FAILED  $folder -> $tag (see $log)" | tee -a "$SUMMARY"
  fi
done
echo "BUILD FINISHED: $(date)" | tee -a "$SUMMARY"
