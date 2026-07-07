#!/usr/bin/env bash
# Build all x86 basej-dnaqc containers locally, one log per image + a summary.
# Arm variants (ubuntu-arm, seqkit-arm, samtools-arm) are intentionally skipped.
set -u

BASE="/home/ubuntu/data/nf-bioskryb-utils/pipelines/basej-dnaqc/container"
LOGDIR="$BASE/_build_logs"
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/SUMMARY.log"
: > "$SUMMARY"

# folder:image_tag  (x86 only)
BUILDS=(
  "ubuntu:basejumper_ubuntu_24.04_stable"
  "seqkit:basejumper_seqkit-2.13.0"
  "samtools:basejumper_samtools-1.23.1"
  "bwamem2_samtools:basejumper_bwamem2_samtools_2.2.1"
  "picard_addorreplacegroups:basejumper_picard_addorreplacegroups_3.1.2"
  "ginkgo:basejumper_ginkgo_0.3.1"
  "preseq_bam2mr:basejumper_preseq_bam2mr_0.1"
  "ginko_parser:basejumper_ginko_parser_0.2.1"
  "custom_parabricks-metrics:basejumper_custom_parabricks-metrics_1.0.3"
  "custom_r_qcplots:basejumper_custom_r_qcplots_0.2.1"
)

echo "BUILD STARTED: $(date)" | tee -a "$SUMMARY"
for entry in "${BUILDS[@]}"; do
  folder="${entry%%:*}"
  tag="${entry##*:}"
  log="$LOGDIR/${folder}.log"
  echo "[$(date +%H:%M:%S)] BUILDING $folder -> $tag (log: $log)" | tee -a "$SUMMARY"
  if docker build -t "$tag" "$BASE/$folder" > "$log" 2>&1; then
    echo "[$(date +%H:%M:%S)] OK      $folder -> $tag" | tee -a "$SUMMARY"
  else
    echo "[$(date +%H:%M:%S)] FAILED  $folder -> $tag (see $log)" | tee -a "$SUMMARY"
  fi
done
echo "BUILD FINISHED: $(date)" | tee -a "$SUMMARY"
