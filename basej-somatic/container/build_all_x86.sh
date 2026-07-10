#!/usr/bin/env bash
# Build all x86 basej-somatic custom containers locally, one log per image + a summary.
# The open-source run starts from a VCF, so the proprietary Sentieon/DNAscope image
# is NOT built here (see README.md).
set -u

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="$BASE/_build_logs"
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/SUMMARY.log"
: > "$SUMMARY"

# folder:image_tag  (x86 only)
BUILDS=(
  "custom_snp_somatic_filter_sequoia_feb2026:basejumper_custom_snp_somatic_filter_sequoia_feb2026"
  "variantannotation_0.4:basejumper_variantannotation_0.4"
  "utility-lite-0.2.0:basejumper_utility-lite-0.2.0"
  "r_ado_1.2:basejumper_r_ado_1.2"
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
