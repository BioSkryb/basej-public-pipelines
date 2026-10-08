nextflow.enable.dsl=2

// ============================================================================
// BASEJ-WGSQC: Rust WGS QC (rewrites.bio rewrite of the basej-wgs FASTQ path)
// ============================================================================
// End-to-end FASTQ -> wgsqc_summary Parquet, mirroring the concept proven in
// basej-dnaqc but targeting the WGS QC metric set (full-genome coverage) that
// basej-wgs collects with the Sentieon driver:
//
//   MERGE_MULTILANE_FASTQ (opt) -> SEQKIT_SAMPLE (opt subsample)
//     -> ALIGN_DEDUP_QC (minibwa-rs | dupblaster | bam-insort + bskryb-qc metrics,
//        one fused task; BAM is optionally published after local QC; WgsMetrics coverage)
//     -> WGSQC_METRICS_TO_PARQUET (wgsqc_summary Parquet, identical Athena schema)
//
// Unlike basej-dnaqc there is NO fastp step (the basej-wgs FASTQ path does not
// trim), and coverage is reported over the full-genome wgs intervals with a
// supplied GENOME_TERRITORY (--genome-territory).
//
// Alignment uses minibwa (BWA-MEM lineage, Rust): NEON on ARM/Graviton by default,
// with an SSE x86 image available for controlled parity runs. Reproduces run
// 1OSLCIDVxdn8WE (BioSkryb/Development, basej-wgs).
// ============================================================================

// ============================================================================
// PROCESS: MERGE_MULTILANE_FASTQ — cat multi-lane R1s / R2s for one biosample
//   Only runs when read1/read2 contain pipe-delimited ("|") multi-lane paths.
// ============================================================================
process MERGE_MULTILANE_FASTQ {
    tag "${sample_name}"

    // Lanes are staged under fixed names and concatenated in CSV order, so merged R1/R2
    // stay pair-synchronised regardless of the original file names (the previous
    // `cat *R1*` glob depended on "R1"/"R2" appearing in the lane file names).
    input:
    tuple val(sample_name), path(r1_files, stageAs: 'lane_R1_*.fastq.gz'), path(r2_files, stageAs: 'lane_R2_*.fastq.gz')

    output:
    tuple val(sample_name), path("${sample_name}_merged_R{1,2}.fastq.gz"), emit: reads

    script:
    def r1_list = (r1_files instanceof List ? r1_files : [r1_files]).join(' ')
    def r2_list = (r2_files instanceof List ? r2_files : [r2_files]).join(' ')
    """
    set -euo pipefail
    cat ${r1_list} > ${sample_name}_merged_R1.fastq.gz & p1=\$!
    cat ${r2_list} > ${sample_name}_merged_R2.fastq.gz & p2=\$!
    wait \$p1
    wait \$p2
    """
}

// ============================================================================
// PROCESS: SEQKIT_SAMPLE — optional subsample + raw read count (total_reads)
//   WGS QC usually runs on full data (skip_subsampling=true). When enabled, caps
//   to max_total_reads with seqkit proportion sampling (same logic as basej-wgs).
//   Always emits read_counts.txt: line1=TOTAL_READS(raw), line2=FINAL_READS.
// ============================================================================
process SEQKIT_SAMPLE {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(max_total_reads)
    val(seqkit_sample_seed)
    val(skip_subsampling)

    output:
    tuple val(sample_name), path("${sample_name}_subsampled_R*.fastq.gz"), emit: reads
    tuple val(sample_name), path("${sample_name}_read_counts.txt"), emit: read_counts

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    set -euo pipefail
    ( zcat '${r1}' | wc -l | awk '{print int(\$1/4)}' > read1.txt ) &
    ( zcat '${r2}' | wc -l | awk '{print int(\$1/4)}' > read2.txt ) &
    wait
    R1=\$(cat read1.txt); R2=\$(cat read2.txt)
    export TOTAL_READS=\$((R1 + R2))

    if [ "${skip_subsampling}" = "true" ] || [ "\$TOTAL_READS" -le "${max_total_reads}" ]; then
        cp '${r1}' '${sample_name}_subsampled_R1.fastq.gz'
        cp '${r2}' '${sample_name}_subsampled_R2.fastq.gz'
        export FINAL_READS=\$TOTAL_READS
    else
        PROPORTION=\$(awk -v t="${max_total_reads}" -v tot="\$TOTAL_READS" 'BEGIN { printf "%.18f", t/tot }')
        seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j ${task.cpus} -o '${sample_name}_subsampled_R1.fastq.gz' '${r1}'
        seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j ${task.cpus} -o '${sample_name}_subsampled_R2.fastq.gz' '${r2}'
        R1_SUB=\$(zcat '${sample_name}_subsampled_R1.fastq.gz' | wc -l | awk '{print int(\$1/4)}')
        R2_SUB=\$(zcat '${sample_name}_subsampled_R2.fastq.gz' | wc -l | awk '{print int(\$1/4)}')
        export FINAL_READS=\$((R1_SUB + R2_SUB))
    fi
    echo "\$TOTAL_READS" > '${sample_name}_read_counts.txt'
    echo "\$FINAL_READS" >> '${sample_name}_read_counts.txt'
    """
}

// ============================================================================
// PROCESS: ALIGN_DEDUP_QC — FUSED align + mark-duplicates + coordinate sort + QC
//   minibwa-rs (BWA-MEM lineage) -> dupblaster (streaming dedup) -> bam-insort
//   (in-RAM coordinate sort + BAI), then bskryb-qc metrics on the LOCAL BAM in the
//   same task.
//
//   Fusing ALIGN_DEDUP and RUSTQC_METRICS avoids staging the full-depth BAM into a
//   second metrics task and lets bskryb-qc use all task CPUs. By default the BAM is
//   deleted after QC; publish_bam=true retains it as an optional reusable output.
//
//   Scoring -A1 -B4 -O6 -E1 matches BWA-MEM (Sentieon) for concordance. The fused task
//   bundles minibwa-rs + dupblaster + bam-insort + samtools + bskryb-qc in one image.
// ============================================================================
process ALIGN_DEDUP_QC {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)
    tuple path(idx_mbw), path(idx_l2b)
    val(platform)
    path(base_intervals)
    // Own subdir: in exome mode this is the same file as base_intervals (name collision).
    path(cov_intervals, stageAs: 'cov_intervals/*')
    tuple path(reference), path(reference_fai)
    path(gc_windows)
    val(publish_bam)

    output:
    tuple val(sample_name), path("${sample_name}.rustqc.json"), emit: metrics
    tuple val(sample_name), path("${sample_name}.dupblaster.tsv"), emit: stats
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam, optional: true

    script:
    def prefix = idx_mbw.baseName
    // minibwa (the long pole) gets all cpus. bam-insort buffers the incoming stream
    // during alignment (its BGZF-decode threads idle-wait on the pipe), then does its
    // parallel in-RAM sort + multithreaded BGZF write + BAI after alignment finishes,
    // so both phases can use all cpus without meaningful overlap/oversubscription.
    def threads = Math.max(1, task.cpus as int)
    // bam-insort record budget. Kept CONSERVATIVE (~25% of task memory): the parallel
    // region-bucket sort gets its speedup from parallelism, NOT from holding everything
    // in RAM, so a smaller budget just spills more (cheap sequential I/O) and stays well
    // clear of the cgroup limit. Headroom is needed for the concurrent minibwa (~11 GB)
    // + dupblaster (~7 GB) upstream in the pipe, phase-2 parallel bucket loads, and
    // allocator overhead (real RSS ran ~1.5x the counted budget, which OOM-killed a
    // 55 GB budget on the 110 GB box). 0.25 * 110 GB => ~27 GB counted, ~40 GB real.
    def insort_mem_gb = Math.max(4, ((task.memory.toGiga() as long) * 0.25) as int)
    // Many small genomic buckets keep bam-insort's phase-2 memory bounded: each worker
    // loads ONE whole bucket to sort it, so peak phase-2 RAM ~= threads * (data/buckets).
    // With ~245 GB of records (minibwa -a emits secondaries -> ~2 records/template),
    // 4096 buckets => ~60 MB/bucket => ~4 GB across 62 workers. (512 buckets OOM-killed
    // phase 2: ~480 MB/bucket * 62 ~= 30 GB + allocator bloat.)
    def insort_buckets = 4096
    // WGS-style coverage block only in wgs mode; exome target coverage comes from
    // Picard CollectHsMetrics (bskryb-qc omits the coverage block without a territory).
    def territory = (params.mode != 'exome' && params.genome_territory) ? "--genome-territory ${params.genome_territory}" : ""
    // Exome targets are ~150 bp, so admit reads that OVERLAP a target (Sentieon/Picard
    // --interval semantics); start-in-target dropped ~40% of on-target reads.
    def overlap_arg = params.mode == 'exome' ? "--interval-overlap" : ""
    def ref_arg = reference.name != 'NO_FILE' ? "--reference ${reference}" : ""
    def gcwin_arg = gc_windows.name != 'NO_FILE.gcwin' ? "--gc-windows ${gc_windows}" : ""
    """
    set -euo pipefail
    export TMPDIR=\$PWD
    # Cap glibc per-thread malloc arenas: with 62+ worker threads churning ~245 GB of
    # small record allocations, the default (8 * ncpu) arenas retain freed memory and
    # bloat RSS. MALLOC_ARENA_MAX=2 keeps RSS close to live memory.
    export MALLOC_ARENA_MAX=2
    # minibwa-rs (align) | dupblaster (streaming dedup) | bam-insort (parallel region-
    # bucket coordinate sort + BAI). bam-insort replaces `samtools sort | samtools index`:
    # it partitions records into genomic buckets, sorts+compresses them in parallel, and
    # concatenates in genomic order — no serial merge. Output validated identical to
    # samtools sort/index. Memory-bounded via --max-mem-gb (spills) + --buckets.
    minibwa-rs map -a -x sr -A1 -B4 -O6 -E1 -Y -t ${threads} \\
        -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:${platform}" \\
        ${prefix} ${reads[0]} ${reads[1]} \\
      | dupblaster --single-end-strategy picard-exact --tmp-dir \$PWD \\
            --stats ${sample_name}.dupblaster.tsv --sample ${sample_name} -o - \\
      | bam-insort -t ${threads} --max-mem-gb ${insort_mem_gb} --buckets ${insort_buckets} --tmp \$PWD -o ${sample_name}.bam

    # Metrics on the LOCAL sorted BAM (no S3 round-trip), region-sharded across all cpus.
    # preseq/Lorenz skipped: low-pass single-cell metrics, meaningless on high-pass WGS.
    bskryb-qc \\
        --bam ${sample_name}.bam --intervals ${base_intervals} \\
        --coverage-intervals ${cov_intervals} --sample ${sample_name} \\
        --threads ${task.cpus} --no-preseq --no-lorenz \\
        ${territory} ${ref_arg} ${gcwin_arg} ${overlap_arg} \\
        --out ${sample_name}.rustqc.json

    if [ "${publish_bam}" = "true" ]; then
        # Optional reusable alignment output. Fail loudly if bam-insort did not produce
        # both members; optional tuple semantics must not hide a partial BAM result.
        test -s ${sample_name}.bam
        test -s ${sample_name}.bam.bai
    else
        rm -f ${sample_name}.bam ${sample_name}.bam.bai
    fi
    """
}

// ============================================================================
// PROCESS: SAMTOOLS_SUBSAMPLE_CRAM — Ultima: optional CRAM subsample to max_total_reads
//   Only runs when skip_subsampling=false (WGS QC normally runs on full data).
//   Counts PRIMARY records; samtools -s keeps whole templates by QNAME hash.
//   Writes the shared read_counts.txt format (line1 TOTAL_READS, line2 FINAL_READS).
// ============================================================================
process SAMTOOLS_SUBSAMPLE_CRAM {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(cram), path(crai), val(target_reads)
    val(samtools_seed)
    tuple path(ref_fasta), path(ref_fai)

    output:
    tuple val(sample_name), path("${sample_name}.sub.cram"), path("${sample_name}.sub.cram.crai"), emit: cram
    tuple val(sample_name), path("${sample_name}_read_counts.txt"), emit: read_counts

    script:
    """
    set -euo pipefail
    count_primary() { samtools view -c -F 0x900 -@ ${task.cpus} --reference '${ref_fasta}' "\$1"; }

    TOTAL_READS=\$(count_primary '${cram}')
    if [ "\$TOTAL_READS" -le "${target_reads}" ]; then
        ln -s '${cram}' '${sample_name}.sub.cram'
        ln -s '${crai}' '${sample_name}.sub.cram.crai'
        FINAL_READS=\$TOTAL_READS
    else
        FRAC=\$(awk -v t="${target_reads}" -v tot="\$TOTAL_READS" 'BEGIN { printf "%.8f", t/tot }' | sed 's/^0//')
        SEED_INT=\$(echo '${samtools_seed}' | cut -d'.' -f1)
        samtools view -s "\${SEED_INT}\${FRAC}" -@ ${task.cpus} --reference '${ref_fasta}' \\
            -C -o '${sample_name}.sub.cram' '${cram}'
        samtools index '${sample_name}.sub.cram' '${sample_name}.sub.cram.crai'
        FINAL_READS=\$(count_primary '${sample_name}.sub.cram')
    fi
    echo "\$TOTAL_READS" > '${sample_name}_read_counts.txt'
    echo "\$FINAL_READS" >> '${sample_name}_read_counts.txt'
    """
}

// ============================================================================
// PROCESS: CRAM_QC — Ultima: QC metrics on a pre-aligned, vendor-deduplicated CRAM
//   Metrics only, no realignment and no re-dedup (same as basej-wgs's
//   SENTIEON/PICARD_METRICS_CRAM): the 0x400 flags from the Ultima demux are
//   trusted. bskryb-qc reads BAM, so the CRAM is decoded to a LOCAL BAM first
//   (already coordinate-sorted; no sort needed), then deleted. samtools flagstat
//   supplies the duplication rate and the primary read count.
//   Single-end reads land in bskryb-qc's UNPAIRED category (>= 0.1.4).
// ============================================================================
process CRAM_QC {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(cram), path(crai)
    path(base_intervals)
    // Own subdir: in exome mode this is the same file as base_intervals (name collision).
    path(cov_intervals, stageAs: 'cov_intervals/*')
    tuple path(reference), path(reference_fai)
    path(gc_windows)
    // Staged in a subdir: same file name as `reference` when run_gcbias=true.
    tuple path(cram_ref, stageAs: 'cram_ref/*'), path(cram_ref_fai, stageAs: 'cram_ref/*')

    output:
    tuple val(sample_name), path("${sample_name}.rustqc.json"), emit: metrics
    tuple val(sample_name), path("${sample_name}.flagstat.txt"), emit: stats
    tuple val(sample_name), path("${sample_name}_cram_read_counts.txt"), emit: read_counts

    script:
    // WGS-style coverage block only in wgs mode; exome target coverage comes from
    // Picard CollectHsMetrics (bskryb-qc omits the coverage block without a territory).
    def territory = (params.mode != 'exome' && params.genome_territory) ? "--genome-territory ${params.genome_territory}" : ""
    // Exome targets are ~150 bp, so admit reads that OVERLAP a target (Sentieon/Picard
    // --interval semantics); start-in-target dropped ~40% of on-target reads.
    def overlap_arg = params.mode == 'exome' ? "--interval-overlap" : ""
    def ref_arg = reference.name != 'NO_FILE' ? "--reference ${reference}" : ""
    def gcwin_arg = gc_windows.name != 'NO_FILE.gcwin' ? "--gc-windows ${gc_windows}" : ""
    """
    set -euo pipefail
    export TMPDIR=\$PWD

    # flagstat has no --reference flag; pass the CRAM reference as an input-format option
    samtools flagstat -@ ${task.cpus} --input-fmt-option reference='${cram_ref}' '${cram}' > ${sample_name}.flagstat.txt

    # Primary read count (QC-pass) -> read_counts.txt (TOTAL == FINAL: no subsampling here)
    PRIMARY=\$(awk '\$4 == "primary" && NF == 4 { print \$1; exit }' ${sample_name}.flagstat.txt)
    echo "\$PRIMARY" > ${sample_name}_cram_read_counts.txt
    echo "\$PRIMARY" >> ${sample_name}_cram_read_counts.txt

    samtools view -b -@ ${task.cpus} --reference '${cram_ref}' -o ${sample_name}.bam '${cram}'
    samtools index -@ ${task.cpus} ${sample_name}.bam

    bskryb-qc \\
        --bam ${sample_name}.bam --intervals ${base_intervals} \\
        --coverage-intervals ${cov_intervals} --sample ${sample_name} \\
        --threads ${task.cpus} --no-preseq --no-lorenz \\
        ${territory} ${ref_arg} ${gcwin_arg} ${overlap_arg} \\
        --out ${sample_name}.rustqc.json

    rm -f ${sample_name}.bam ${sample_name}.bam.bai
    """
}

// ============================================================================
// PROCESS: PICARD_COLLECTHSMETRICS — exome mode only: hybrid-selection metrics
//   Same command as basej-wgs (bait == target == the panel interval_list). Reads the
//   ALIGN_DEDUP_QC BAM (FASTQ path) or the pre-aligned CRAM (Ultima path); Picard
//   decodes CRAM with REFERENCE_SEQUENCE. bskryb-qc has no HsMetrics module, so this
//   is the source of every target-capture column and of the spec_wes gates.
//   The reference is staged as fasta + fai + dict (Picard needs the .dict alongside).
// ============================================================================
process PICARD_COLLECTHSMETRICS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(aln), path(aln_index)
    tuple path(ref_fasta, stageAs: 'ref/*'), path(ref_fai, stageAs: 'ref/*'), path(ref_dict, stageAs: 'ref/*')
    path(target_intervals)

    output:
    tuple val(sample_name), path("${sample_name}.hsmetrics.txt"), emit: metrics

    script:
    def avail_mem = Math.max(2, (task.memory.toGiga() as int) - 1)
    """
    picard -Xmx${avail_mem}g CollectHsMetrics \\
        --INPUT ${aln} \\
        --OUTPUT ${sample_name}.hsmetrics.txt \\
        --BAIT_INTERVALS ${target_intervals} \\
        --TARGET_INTERVALS ${target_intervals} \\
        --REFERENCE_SEQUENCE ${ref_fasta}
    """
}

// ============================================================================
// PROCESS: GCBIAS_INDEX — one-time windowsByGc precompute (reused by all samples)
// ============================================================================
process GCBIAS_INDEX {
    tag "gc_index"

    input:
    tuple path(reference), path(reference_fai)
    path intervals

    output:
    path("gcwin_w100.txt"), emit: index

    script:
    """
    bskryb-qc gc-index --reference ${reference} --intervals ${intervals} --window 100 --out gcwin_w100.txt
    """
}

// ============================================================================
// PROCESS: RUSTQC_METRICS — folded into ALIGN_DEDUP_QC (runs on the local sorted BAM
//   in the same task; see above). Kept as a comment to document the change.
// ============================================================================

// ============================================================================
// PROCESS: WGSQC_METRICS_TO_PARQUET — merge rustqc + dupblaster + read counts
//          into the canonical wgsqc_summary Parquet and aggregation rows.
// ============================================================================
process WGSQC_METRICS_TO_PARQUET {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(rustqc_json), path(dupblaster_stats), path(read_counts), path(hsmetrics)
    val(mode)
    val(genome)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(user)

    output:
    // wgsqc_* in wgs mode, wesqc_* in exome mode (same table schema, as in basej-wgs)
    path("*qc_summary/workspace=*/workflow_id=*/biosample=*/output.parquet"), emit: parquet
    path("${sample_name}_*qc_metrics.tsv"), emit: metrics_row
    path("${sample_name}_status.csv"), emit: status_row

    script:
    // FASTQ path stages dupblaster stats; the Ultima CRAM path stages a samtools
    // flagstat of the vendor-deduplicated CRAM instead (no re-dedup, as in basej-wgs).
    def dedup_arg = dupblaster_stats.name.endsWith('.flagstat.txt')
        ? "--flagstat ${dupblaster_stats}"
        : "--dupblaster-stats ${dupblaster_stats}"
    def hs_arg = hsmetrics.name != 'NO_HSMETRICS' ? "--hsmetrics ${hsmetrics}" : ""
    """
    wgsqc_metrics_to_parquet.py \\
        --json ${rustqc_json} \\
        ${dedup_arg} ${hs_arg} \\
        --read-counts ${read_counts} \\
        --sample ${sample_name} \\
        --mode ${mode} --genome ${genome} \\
        --dataset-id ${dataset_id} --workspace ${workspace} --workflow-id ${workflow_id} \\
        --pipeline-version ${pipeline_version} --user ${user}
    """
}

// ============================================================================
// PROCESS: AGGREGATE_WGSQC_OUTPUTS — one canonical TSV and platform status index
// ============================================================================
process AGGREGATE_WGSQC_OUTPUTS {
    tag "aggregate_wgsqc_outputs"

    input:
    path(metrics_rows)
    path(status_rows)
    val(mode)

    output:
    path("${mode == 'exome' ? 'wes' : 'wgs'}qc_all_metrics.tsv"), emit: summary_tsv
    path("per_biosample_status.csv"), emit: per_biosample_status

    script:
    def prefix = mode == 'exome' ? 'wes' : 'wgs'
    """
    python3 - <<'PYEOF'
import glob
import pandas as pd

metric_paths = sorted(glob.glob("*_${prefix}qc_metrics.tsv"))
status_paths = sorted(glob.glob("*_status.csv"))
if not metric_paths or not status_paths:
    raise SystemExit("ERROR: no WGSQC metric/status rows were staged for aggregation")

metrics = pd.concat(
    [pd.read_csv(path, sep="\\t", keep_default_na=True) for path in metric_paths],
    ignore_index=True,
).sort_values("biosample")
statuses = pd.concat(
    [pd.read_csv(path, keep_default_na=True) for path in status_paths],
    ignore_index=True,
).sort_values("biosampleName")

if metrics["biosample"].duplicated().any() or statuses["biosampleName"].duplicated().any():
    raise SystemExit("ERROR: duplicate biosample rows in WGSQC aggregation inputs")
if set(metrics["biosample"]) != set(statuses["biosampleName"]):
    raise SystemExit("ERROR: WGSQC metric and status biosample sets do not match")

metrics.to_csv("${prefix}qc_all_metrics.tsv", sep="\\t", index=False, na_rep="NA")
statuses.to_csv("per_biosample_status.csv", index=False)
PYEOF
    """
}

// ============================================================================
// PROCESS: WGS_QC_PLOTS — composition plot + consensus-score tables (shared R)
// ----------------------------------------------------------------------------
// Same R scripts and outputs as basej-wgs's WGS_QC_PLOTS (custom_r_qcplots image):
//   wgs mode  : wgs_qc_plot.R -> qc_wgs.pdf, WGS-QC_composition_mqc.jpg, WGS-QC_*
//   exome mode: wes_qc_plot.R -> qc_wes.pdf, WES-QC_composition_mqc.jpg, WES-QC_*
// The aggregated *qc_all_metrics.tsv already uses the column names these scripts
// read (biosample, total_reads, pct_duplication, pct_chimeras, pct_*x / HS columns),
// so it is passed straight through. The Python verdict (qc_scoring.py) remains the
// authoritative qc_status in Parquet / per_biosample_status; the R CompositeScore is
// checked against it in qc_score_concordance.tsv.
// ============================================================================
process WGS_QC_PLOTS {
    tag "wgs_qc_plots"

    input:
    path(all_metrics_tsv)
    val(mode)

    output:
    path("qc_${mode == 'exome' ? 'wes' : 'wgs'}.pdf"), emit: plots_pdf
    path("*-QC_composition_mqc.jpg"), emit: plots_jpg
    path("*-QC_ConsensusScores.txt"), emit: scores
    path("*-QC_ConsensusScores_SummaryTable_mqc.txt"), emit: scores_summary
    path("*-QC_QCBand_SummaryTable_mqc.txt"), optional: true, emit: scores_bands
    path("qc_score_concordance.tsv"), emit: concordance

    script:
    def prefix = mode == 'exome' ? 'WES' : 'WGS'
    def r_script = mode == 'exome' ? 'wes_qc_plot.R' : 'wgs_qc_plot.R'
    """
    set -euo pipefail

    # Physical copy: R reads it, and the staged symlink must not be rewritten.
    cp ${all_metrics_tsv} input_metrics.tsv

    Rscript /usr/local/bin/${r_script} --metrics_file input_metrics.tsv

    compare_qc_scores.py \\
        --summary input_metrics.tsv \\
        --scores ${prefix}-QC_ConsensusScores.txt \\
        --out qc_score_concordance.tsv

    # MultiQC custom-content header so the score distribution renders as a table.
    tmp=\$(mktemp)
    printf '# id: "qc_score_distribution"\\n# section_name: "Total Usable Cells"\\n# description: "Distribution of samples across composite QC score categories. Scores range from 0 (poor quality) to 5 (high quality)."\\n# plot_type: "table"\\n# pconfig:\\n#   id: "qc_score_dist_table"\\n#   title: "Total Usable Cells"\\n' \\
        | cat - ${prefix}-QC_ConsensusScores_SummaryTable_mqc.txt > "\$tmp"
    mv "\$tmp" ${prefix}-QC_ConsensusScores_SummaryTable_mqc.txt
    """
}

// ============================================================================
// PROCESS: MULTIQC — user-facing WGS QC summary report
// ============================================================================
process MULTIQC {
    tag "multiqc"

    input:
    path(metrics_tsv)
    path(status_csv)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(nextflow_version)
    val(architecture)
    val(fused_container)
    val(seqkit_container)
    path(qc_plot_files)
    path(logo)
    val(mode)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data
    path("tool_mqc_versions.yml"), emit: versions

    script:
    def exome = mode == 'exome'
    def label = exome ? 'WES' : 'WGS'
    def picard_versions = exome
        ? '  picard-hsmetrics: "3.0.0"\\n  picard-container: "quay.io/biocontainers/picard:3.0.0--hdfd78af_0"\\n'
        : ''
    """
    set -euo pipefail

    MULTIQC_VERSION=\$(multiqc --version 2>&1 | head -n 1 | sed -E 's/^multiqc, version[[:space:]]+//' | tr -d '\r')
    test -n "\$MULTIQC_VERSION"

    # Legacy BioSkryb/MultiQC convention: this filename and grouped YAML shape are
    # auto-discovered and rendered in the standard Software Versions section.
    cat > tool_mqc_versions.yml << EOF
basej-wgsqc:
  pipeline: "${pipeline_version}"
  nextflow: "${nextflow_version}"
  minibwa-rs: "0.6-r416"
  dupblaster: "0.1.1"
  bam-insort: "0.1.0"
  bskryb-qc: "0.1.5"
  samtools: "1.16.1"
  seqkit: "2.13.0"
  multiqc: "\$MULTIQC_VERSION"
  execution-architecture: "${architecture}"
  fused-container: "${fused_container}"
  seqkit-container: "${seqkit_container}"
  multiqc-container: "quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
EOF
    printf '${picard_versions}' >> tool_mqc_versions.yml

    python3 - <<'PYEOF'
import csv

with open("${metrics_tsv}", newline="") as handle:
    metrics = list(csv.DictReader(handle, delimiter="\\t"))
with open("${status_csv}", newline="") as handle:
    statuses = {row["biosampleName"]: row["qc_status"] for row in csv.DictReader(handle)}

if not metrics:
    raise SystemExit("ERROR: no WGSQC metrics rows were staged for MultiQC")
if {row["biosample"] for row in metrics} != set(statuses):
    raise SystemExit("ERROR: WGSQC metrics and status biosample sets do not match")
for row in metrics:
    if row["qc_status"] != statuses[row["biosample"]]:
        raise SystemExit(f"ERROR: QC status mismatch for {row['biosample']}")

exome = "${mode}" == "exome"
label = "${label}"
# Exome columns match basej-wgs's wesqc MultiQC table (target capture instead of WGS coverage).
columns = [
    ("sample_name", "biosample"),
    ("QC_Status", "qc_status"),
    ("QC_Score", "qc_score"),
    ("Total_Reads", "total_reads"),
    ("Align_Total_Reads", "align_total_reads"),
    ("PCT_Target_10x", "pct_target_bases_10x"),
    ("Zero_Cvg_Targets_Pct", "zero_cvg_targets_pct"),
    ("Fold_80_Base_Penalty", "fold_80_base_penalty"),
    ("Mean_Target_Coverage", "mean_target_coverage"),
    ("PCT_Selected_Bases", "pct_selected_bases"),
    ("PCT_Duplication", "pct_duplication"),
    ("Insert_Median", "insert_median"),
    ("PCT_Chimeras", "pct_chimeras"),
] if exome else [
    ("sample_name", "biosample"),
    ("QC_Status", "qc_status"),
    ("QC_Score", "qc_score"),
    ("Total_Reads", "total_reads"),
    ("Align_Total_Reads", "align_total_reads"),
    ("PCT_Duplication", "pct_duplication"),
    ("PCT_Chimeras", "pct_chimeras"),
    ("Mean_Coverage", "mean_coverage"),
    ("PCT_1x", "pct_1x"),
    ("PCT_5x", "pct_5x"),
    ("PCT_10x", "pct_10x"),
    ("PCT_30x", "pct_30x"),
    ("Insert_Median", "insert_median"),
    ("AT_Dropout", "at_dropout"),
    ("GC_Dropout", "gc_dropout"),
]

with open("wgsqc_summary_mqc.txt", "w", newline="") as handle:
    handle.write("# id: 'wgsqc_summary'\\n")
    handle.write("# plot_type: 'table'\\n")
    handle.write("# section_name: 'QC Summary'\\n")
    handle.write("# description: 'Per-sample alignment, "
                 + ("target capture" if exome else "coverage")
                 + ", duplication, and QC status.'\\n")
    handle.write("# pconfig:\\n")
    handle.write("#   id: 'wgsqc_summary_table'\\n")
    handle.write(f"#   title: '{label} QC Summary'\\n")
    # lineterminator="\\n": csv defaults to "\\r\\n", which against the newline="" handle above
    # would emit CRLF data rows under LF comment lines and leave a stray CR on the last
    # column's value for MultiQC to parse.
    writer = csv.DictWriter(handle, fieldnames=[name for name, _ in columns], delimiter="\\t",
                            lineterminator="\\n")
    writer.writeheader()
    for row in sorted(metrics, key=lambda item: item["biosample"]):
        writer.writerow({name: row.get(source, "NA") for name, source in columns})
PYEOF

    cat > multiqc_config.yaml << EOF
custom_logo_title: 'BioSkryb Genomics'
custom_logo: ${logo}
custom_logo_width: 260

title: "basej-${exome ? 'wesqc' : 'wgsqc'}-rs v${pipeline_version}"
report_header_info:
  - Dataset ID: "${dataset_id}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false
skip_generalstats: true

module_order:
  - custom_content

report_section_order:
  wgsqc_summary:
    order: 1000
  WGS-QC_composition:
    order: 900
  WES-QC_composition:
    order: 900
  qc_score_distribution:
    order: 800

custom_data:
  WGS-QC_composition:
    section_name: "QC Composition"
    description: "Per-sample WGS QC metric distributions (coverage breadth, duplication, chimeras) grouped by QC tier."
  WES-QC_composition:
    section_name: "QC Composition"
    description: "Per-sample exome QC metric distributions (target coverage, uniformity, bait efficiency) grouped by QC tier."

table_cond_formatting_rules:
  QC_Status:
    pass:
      - s_eq: "PASS"
    warn:
      - s_eq: "Borderline"
    fail:
      - s_eq: "FAIL"
EOF

    # tool_mqc_versions.yml is staged as an input and auto-discovered by MultiQC.
    multiqc . -n multiqc_report.html -c multiqc_config.yaml --force
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================
workflow {
    main:
    if (!(params.mode in ['wgs', 'exome'])) {
        error "ERROR: params.mode must be 'wgs' or 'exome'; got '${params.mode}'"
    }
    def exome = params.mode == 'exome'
    if (!exome && (!params.genome_territory || params.genome_territory.toString().toLong() <= 0L)) {
        error "basej-wgsqc requires a positive genome_territory for WGS coverage metrics"
    }

    // Exome (WES) mode: AlignmentStat/GcBias/InsertSize AND target coverage all run over the
    // capture panel's target intervals (same as basej-wgs's exome driver), and target-capture
    // metrics come from Picard CollectHsMetrics. --target_intervals overrides the panel lookup.
    def target_intervals = null
    if (exome) {
        def panel = params.genomes[params.genome]?.get(params.exome_panel)
        target_intervals = params.target_intervals ?: panel?.get('wgs_or_target_intervals')
        if (!target_intervals) {
            error "ERROR: mode='exome' needs target intervals: exome_panel '${params.exome_panel}' " +
                  "has no wgs_or_target_intervals for genome '${params.genome}' in genomes.config " +
                  "(or pass --target_intervals <panel.interval_list>)."
        }
        log.info "basej-wgsqc: exome mode, panel '${params.exome_panel}', targets ${target_intervals}"
    }

    def resolved_pipeline_version = params.pipeline_version ?: workflow.manifest.version
    def resolved_nextflow_version = workflow.nextflow.version.toString()
    def fused_container = "basejumper_minibwa-rustqc_0.1.5"
    def seqkit_container = "basejumper_seqkit-2.13.0"

    // Input CSV, one row per biosample. Rows are routed by content (same rule as basej-wgs),
    // so a run may mix both kinds:
    //   FASTQ (Illumina/Element): biosampleName (or biosample), read1, read2
    //         read1/read2 may hold "|"-delimited multi-lane paths.
    //   CRAM  (Ultima, pre-aligned + vendor-deduplicated): biosampleName, cram[, crai]
    ch_rows = Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def name = (row.biosampleName ?: row.biosample)?.trim()
            if (!name) {
                error "ERROR: input_csv row is missing biosampleName: ${row}"
            }
            [name, row]
        }
        .branch { name, row ->
            cram: row.cram?.trim()
            fastq: row.read1?.trim() && row.read2?.trim()
            other: true
        }
    ch_rows.other.subscribe { name, row ->
        error "ERROR: biosample '${name}' needs either read1+read2 or cram in input_csv"
    }

    // Base-metrics set drives AlignmentStat/GcBias/InsertSize; the narrower coverage set
    // drives only the WgsMetrics-style coverage block (see nextflow.config).
    // Exome: both sets are the panel targets.
    ch_base_intervals = Channel.value(file(exome ? target_intervals : params.intervals))
    ch_cov_intervals  = Channel.value(file(exome ? target_intervals
                                                 : (params.coverage_intervals ?: params.wgs_or_target_intervals)))
    ch_minibwa_index = Channel.value([file("${params.minibwa_index}.mbw"), file("${params.minibwa_index}.l2b")])

    // GcBias reference + windows index (enabled by default for Sentieon parity).
    ch_reference = params.run_gcbias
        ? Channel.value([file(params.reference_fasta), file("${params.reference_fasta}.fai")])
        : Channel.value([file("${projectDir}/assets/NO_FILE"), file("${projectDir}/assets/NO_FILE.fai")])
    if (params.run_gcbias) {
        // The precomputed windowsByGc index matches the WGS base intervals only; in exome
        // mode it is rebuilt once over the panel targets.
        if (params.gc_windows_index && !exome) {
            ch_gcwin = Channel.value(file(params.gc_windows_index))
        } else {
            GCBIAS_INDEX(ch_reference, ch_base_intervals)
            ch_gcwin = GCBIAS_INDEX.out.index.first()
        }
    } else {
        ch_gcwin = Channel.value(file("${projectDir}/assets/NO_FILE.gcwin"))
    }

    // ===== FASTQ path: (merge lanes) -> (opt subsample) -> fused align/dedup/sort/QC =====
    ch_fastq = ch_rows.fastq.map { name, row ->
        def r1 = row.read1.tokenize('|').collect { file(it.trim()) }
        def r2 = row.read2.tokenize('|').collect { file(it.trim()) }
        if (r1.size() != r2.size()) {
            error "ERROR: biosample '${name}' has ${r1.size()} read1 lane(s) but ${r2.size()} read2 lane(s)"
        }
        tuple(name, r1, r2)
    }
    ch_fastq_branched = ch_fastq.branch { name, r1, r2 ->
        multilane: r1.size() > 1
        singlelane: true
    }
    MERGE_MULTILANE_FASTQ(ch_fastq_branched.multilane)
    ch_reads = ch_fastq_branched.singlelane
        .map { name, r1, r2 -> tuple(name, [r1[0], r2[0]]) }
        .mix(MERGE_MULTILANE_FASTQ.out.reads)

    SEQKIT_SAMPLE(
        ch_reads.map { s, r -> tuple(s, r, params.max_total_reads) },
        params.seqkit_sample_seed,
        params.skip_subsampling
    )

    // Fused align/dedup/sort/QC. The BAM is kept when it is published (publish_bam, default
    // true) and always in exome mode, where CollectHsMetrics reads it.
    ALIGN_DEDUP_QC(
        SEQKIT_SAMPLE.out.reads, ch_minibwa_index, params.platform ?: 'ILLUMINA',
        ch_base_intervals, ch_cov_intervals, ch_reference, ch_gcwin,
        params.publish_bam || exome
    )

    // ===== CRAM path (Ultima): (opt subsample) -> metrics on the vendor alignment =====
    ch_cram = ch_rows.cram.map { name, row ->
        def cram = file(row.cram.trim())
        def crai = row.crai?.trim() ? file(row.crai.trim()) : file("${row.cram.trim()}.crai")
        tuple(name, cram, crai)
    }
    ch_cram_ref = Channel.value([
        file(params.reference_fasta),
        file("${params.reference_fasta}.fai")
    ])
    if (params.skip_subsampling) {
        ch_cram_for_qc = ch_cram
    } else {
        SAMTOOLS_SUBSAMPLE_CRAM(
            ch_cram.map { s, c, i -> tuple(s, c, i, params.max_total_reads) },
            params.samtools_seed,
            ch_cram_ref
        )
        ch_cram_for_qc = SAMTOOLS_SUBSAMPLE_CRAM.out.cram
    }
    CRAM_QC(
        ch_cram_for_qc, ch_base_intervals, ch_cov_intervals, ch_reference, ch_gcwin, ch_cram_ref
    )
    // Raw/final read counts: from the subsample step when it ran (raw = full CRAM),
    // otherwise the primary count of the CRAM that was measured.
    ch_cram_read_counts = params.skip_subsampling
        ? CRAM_QC.out.read_counts
        : SAMTOOLS_SUBSAMPLE_CRAM.out.read_counts

    // Join per-sample metric inputs by biosample for the Parquet writer.
    // FASTQ rows carry dupblaster stats; CRAM rows carry samtools flagstat.
    ch_merged_base = ALIGN_DEDUP_QC.out.metrics
        .join(ALIGN_DEDUP_QC.out.stats)
        .join(SEQKIT_SAMPLE.out.read_counts)
        .mix(
            CRAM_QC.out.metrics
                .join(CRAM_QC.out.stats)
                .join(ch_cram_read_counts)
        )

    // ===== Exome: Picard CollectHsMetrics on the BAM (FASTQ path) or CRAM (Ultima) =====
    if (exome) {
        def ref_dict = params.reference_fasta.toString().replaceAll(/\.fa(sta)?$/, '.dict')
        PICARD_COLLECTHSMETRICS(
            ALIGN_DEDUP_QC.out.bam.mix(ch_cram_for_qc),
            Channel.value([
                file(params.reference_fasta),
                file("${params.reference_fasta}.fai"),
                file(ref_dict)
            ]),
            Channel.value(file(target_intervals))
        )
        ch_merged = ch_merged_base.join(PICARD_COLLECTHSMETRICS.out.metrics)
    } else {
        def no_hs = file("${projectDir}/assets/NO_HSMETRICS")
        ch_merged = ch_merged_base.map { row -> row + [no_hs] }
    }

    WGSQC_METRICS_TO_PARQUET(
        ch_merged,
        params.mode, params.genome,
        params.dataset_id, params.workspace, params.workflow_id,
        resolved_pipeline_version, params.pipeline_user
    )

    AGGREGATE_WGSQC_OUTPUTS(
        WGSQC_METRICS_TO_PARQUET.out.metrics_row.collect(),
        WGSQC_METRICS_TO_PARQUET.out.status_row.collect(),
        params.mode
    )

    WGS_QC_PLOTS(AGGREGATE_WGSQC_OUTPUTS.out.summary_tsv, params.mode)
    ch_qc_plots = WGS_QC_PLOTS.out.plots_pdf
        .mix(WGS_QC_PLOTS.out.plots_jpg)
        .mix(WGS_QC_PLOTS.out.scores_bands)
        .mix(WGS_QC_PLOTS.out.concordance)

    MULTIQC(
        AGGREGATE_WGSQC_OUTPUTS.out.summary_tsv,
        AGGREGATE_WGSQC_OUTPUTS.out.per_biosample_status,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        resolved_pipeline_version,
        resolved_nextflow_version,
        params.architecture,
        fused_container,
        seqkit_container,
        WGS_QC_PLOTS.out.plots_jpg.mix(WGS_QC_PLOTS.out.scores_summary).collect(),
        file("${projectDir}/assets/bioskryb_logo-tagline.png", checkIfExists: true),
        params.mode
    )

    publish:
    // publish_bam=false in exome mode still keeps the BAM for HsMetrics; don't publish it.
    bam_files = (params.publish_bam ? ALIGN_DEDUP_QC.out.bam : channel.empty())
        .map { sample_name, bam, bai ->
            [biosampleName: sample_name, bam: bam, bai: bai]
        }
    rustqc_metrics = ALIGN_DEDUP_QC.out.metrics
        .mix(CRAM_QC.out.metrics)
        .map { sample_name, metrics ->
            [biosampleName: sample_name, metrics: metrics]
        }
    // FASTQ: dupblaster stats; Ultima CRAM: samtools flagstat of the vendor-deduplicated CRAM
    dedup_metrics = ALIGN_DEDUP_QC.out.stats
        .mix(CRAM_QC.out.stats)
        .map { sample_name, metrics ->
            [biosampleName: sample_name, metrics: metrics]
        }
    wgsqc_summary = WGSQC_METRICS_TO_PARQUET.out.parquet
    wgsqc_summary_tsv = AGGREGATE_WGSQC_OUTPUTS.out.summary_tsv
    per_biosample_status = AGGREGATE_WGSQC_OUTPUTS.out.per_biosample_status
    qc_plots = ch_qc_plots
    // Consensus-score tables: metrics/{wgs,wes}qc_metrics/, same as basej-wgs's wgsqc_scores
    wgsqc_scores = WGS_QC_PLOTS.out.scores.mix(WGS_QC_PLOTS.out.scores_summary)
    software_versions = MULTIQC.out.versions
    multiqc_report = MULTIQC.out.report
}

// ============================================================================
// PLATFORM OUTPUTS
// ============================================================================
output {
    bam_files {
        path "bam/${params.workspace}/dna/tool=minibwa-rs/pipeline=wgsqc"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/bam.csv"
            header true
        }
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "bam",
             tool: "minibwa-rs",
             reference: params.genome
    }

    rustqc_metrics {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/rustqc"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/metrics.csv"
            header true
        }
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "qc_metrics",
             tool: "bskryb-qc",
             reference: params.genome
    }

    dedup_metrics {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/dedup"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/dedup_metrics.csv"
            header true
        }
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "dedup_metrics",
             tool: "dupblaster",
             reference: params.genome
    }

    wgsqc_summary {
        path "tables"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_summary".toString(),
             reference: params.genome
    }

    wgsqc_summary_tsv {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/${params.mode == 'exome' ? 'wes' : 'wgs'}qc_metrics"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_summary_tsv".toString()
    }

    wgsqc_scores {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/${params.mode == 'exome' ? 'wes' : 'wgs'}qc_metrics"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_scores".toString()
    }

    // Composition plot, QC band table and Python/R concordance (basej-wgs qc_plots layout)
    qc_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/qc_plots"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_plots".toString()
    }

    per_biosample_status {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/index"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             artifact: "per_biosample_status"
    }

    software_versions {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/execution_info"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             artifact: "software_versions"
    }

    multiqc_report {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/reports"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             artifact: "multiqc_report"
    }
}
