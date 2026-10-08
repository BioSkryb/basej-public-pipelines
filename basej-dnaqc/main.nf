nextflow.enable.dsl=2

// ============================================================================
// BASEJ-DNAQC: Rust DNA QC (rewrites.bio rewrite of basej-dnaqc)
// ============================================================================
// End-to-end FASTQ -> dnaqc_summary Parquet:
//   SEQKIT_SAMPLE (subsample 2M) -> FASTP_TRIM
//   -> ALIGN_DEDUP (minibwa-rs | dupblaster | samtools sort, one fused pipe)
//   -> RUSTQC_METRICS (single-pass Rust metrics) -> METRICS_TO_PARQUET
//
// Alignment uses minibwa (BWA-MEM lineage, Rust). Runs on Graviton (arm64,
// default) via minibwa's s2n_lite scalar SW path — bit-identical to the x86 SSE
// build; a NEON port is a later perf optimization. All BAM-derived metrics come
// from the single-pass Rust binary. Set architecture='x86' for the SSE images.
//
// CNV (Ginkgo) + MultiQC parity with basej-methylqc-rs: GINKO_NOPUBLISH runs the
// Ginkgo CNV chain (unchanged), GINKGO_BINS_TO_PARQUET writes the bin-level
// cnv_summary Parquet, METRICS_TO_PARQUET is an aggregate collect-based writer, and
// MULTIQC renders the QC report.
// ============================================================================

// Ginkgo CNV workflow (cloned config from basej-dnaqc/modules.nf)
include { GINKO_NOPUBLISH } from './modules.nf'

// ============================================================================
// PROCESS: MERGE_MULTILANE_FASTQ — concatenate multi-lane R1s / R2s per biosample
//   Only runs for rows whose read1/read2 contain pipe-delimited ("|") paths.
//   Lanes are staged under fixed names and concatenated in CSV order (gzip
//   members concatenate into a valid gzip stream), so the merged R1/R2 stay
//   pair-synchronised regardless of the original file names.
// ============================================================================
process MERGE_MULTILANE_FASTQ {
    tag "${sample_name}"

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
// PROCESS: SEQKIT_SAMPLE — subsample to target reads (2M) for QC
// ============================================================================
process SEQKIT_SAMPLE {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(target_reads)
    val(seqkit_sample_seed)

    output:
    tuple val(sample_name), path("${sample_name}_subsampled_R*.fastq.gz"), emit: reads
    tuple val(sample_name), path("${sample_name}_read_counts.txt"), emit: read_counts

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    // R1 and R2 are sampled CONCURRENTLY, splitting the task's cores between them.
    // This was the dominant cost of the pipeline: on run 1DaSflsTOdCaHd the two
    // sequential `seqkit sample` calls plus their sequential re-counts held the task
    // at ~1.1 busy cores for 2.5-4.5 h per sample. Selection is per-file and seeded
    // (-s), and seqkit draws per record in a single streaming pass, so running the
    // mates in parallel does not change which reads are picked -- pairs stay in sync.
    def sample_threads = Math.max(1, (task.cpus as int).intdiv(2))
    """
    set -euo pipefail

    # zcat|wc -l over a full FASTQ is the expensive part; always run the pair in
    # parallel and propagate failures (bare `wait` returns 0 and would hide them,
    # silently yielding a wrong total_reads).
    count_reads() { zcat "\$1" | wc -l | awk '{print int(\$1/4)}'; }

    count_reads '${r1}' > read1.txt & p1=\$!
    count_reads '${r2}' > read2.txt & p2=\$!
    wait \$p1
    wait \$p2

    R1=\$(cat read1.txt); R2=\$(cat read2.txt)
    export TOTAL_READS=\$((R1 + R2))
    if [ "\$TOTAL_READS" -le "${target_reads}" ]; then
        cp '${r1}' '${sample_name}_subsampled_R1.fastq.gz'
        cp '${r2}' '${sample_name}_subsampled_R2.fastq.gz'
        export FINAL_READS=\$TOTAL_READS
    else
        PROPORTION=\$(awk -v t="${target_reads}" -v tot="\$TOTAL_READS" 'BEGIN { printf "%.18f", t/tot }')
        seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j ${sample_threads} -o '${sample_name}_subsampled_R1.fastq.gz' '${r1}' & s1=\$!
        seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j ${sample_threads} -o '${sample_name}_subsampled_R2.fastq.gz' '${r2}' & s2=\$!
        wait \$s1
        wait \$s2

        count_reads '${sample_name}_subsampled_R1.fastq.gz' > sub1.txt & q1=\$!
        count_reads '${sample_name}_subsampled_R2.fastq.gz' > sub2.txt & q2=\$!
        wait \$q1
        wait \$q2

        export FINAL_READS=\$((\$(cat sub1.txt) + \$(cat sub2.txt)))
    fi
    echo "\$TOTAL_READS" > '${sample_name}_read_counts.txt'
    echo "\$FINAL_READS" >> '${sample_name}_read_counts.txt'
    """
}

// ============================================================================
// PROCESS: SAMTOOLS_SUBSAMPLE_CRAM — Ultima: subsample a pre-aligned CRAM to n_reads
//   Same role as SEQKIT_SAMPLE for the FASTQ path (and same read_counts.txt format:
//   line1 = TOTAL_READS, line2 = FINAL_READS). Counts PRIMARY records only
//   (-F 0x900), so a supplementary alignment is not counted as an extra read.
//   samtools -s keeps/drops whole templates by QNAME hash, so supplementaries
//   follow their primary. Seeded and deterministic.
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
        # samtools -s takes SEED.FRACTION (integer seed, fractional keep-probability)
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
// PROCESS: CRAM_DEDUP — Ultima: decode the (subsampled) pre-aligned CRAM and re-mark
//   duplicates, without realignment (matches basej-dnaqc's CRAM path, which keeps
//   the vendor alignment and re-deduplicates after subsampling).
//
//   Re-deduplication is required: the vendor 0x400 flags were computed on the full-
//   depth library, so after subsampling most "duplicates" no longer have their
//   original in the set (measured: 7.4% flagged vs 0.18% real on a 500K-read Ultima
//   subsample). samtools collate -> query-grouped stream -> dupblaster (strand-aware
//   single-end key, the Picard fragSort equivalent) -> samtools sort. Output shape is
//   identical to ALIGN_DEDUP (BAM+BAI, dupblaster stats TSV).
// ============================================================================
process CRAM_DEDUP {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(cram), path(crai)
    tuple path(ref_fasta), path(ref_fai)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("${sample_name}.dupblaster.tsv"), emit: stats

    script:
    def sort_threads = Math.max(1, (task.cpus as int).intdiv(2))
    def sort_mem_mb  = Math.max(768, ((task.memory.toMega() as long) * 0.4 / sort_threads) as int)
    """
    set -euo pipefail
    export TMPDIR=\$PWD
    samtools collate -O -u -@ ${task.cpus} --reference '${ref_fasta}' '${cram}' collate_tmp \\
      | dupblaster --single-end-strategy strand-aware --tmp-dir \$PWD \\
            --stats ${sample_name}.dupblaster.tsv --sample ${sample_name} -o - \\
      | samtools sort -@ ${sort_threads} -m ${sort_mem_mb}M -T ${sample_name}.sort -o ${sample_name}.bam -
    samtools index -@ ${task.cpus} ${sample_name}.bam
    """
}

// ============================================================================
// PROCESS: FASTP_TRIM — adapter trimming + read-level QC JSON (fastp_* fields)
// ============================================================================
process FASTP_TRIM {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*_trim.fastq.gz"), emit: reads
    tuple val(sample_name), path("${sample_name}_fastp.json"), emit: json

    script:
    """
    fastp --thread ${task.cpus} \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${sample_name}_R1_trim.fastq.gz --out2 ${sample_name}_R2_trim.fastq.gz \\
        --json ${sample_name}_fastp.json --html ${sample_name}_fastp.html \\
        --detect_adapter_for_pe \\
        2> ${sample_name}_fastp.log
    """
}

// ============================================================================
// PROCESS: ALIGN_DEDUP — single-pass align + mark-duplicates + coordinate sort
//   minibwa-rs (BWA-MEM lineage, x86) -> dupblaster (streaming dedup) -> samtools sort
//
//   Everything runs in one process over a single streaming pipe, so the BAM is
//   materialized exactly once (the final coordinate-sorted, indexed BAM). No
//   intermediate unsorted/aligned BAM is staged between Nextflow tasks.
//
//   - minibwa-rs `-t` threads emit query-grouped SAM/BAM (R1/R2 contiguous), which
//     is exactly what dupblaster consumes with no coordinate sort.
//   - dupblaster (Fulcrum Genomics, MIT; Rust successor to samblaster / Picard
//     MarkDuplicates) marks dups in one pass and writes uncompressed BAM (level 0)
//     straight into the sort, avoiding a recompress round-trip. `picard-exact`
//     gives 100% Picard orphan concordance (order-independent; we sort after).
//     NOTE: dupblaster does not detect optical duplicates (documented drift).
//   - samtools sort `-@` threads + `-m` per-thread memory for the coordinate sort.
//
//   Scoring -A1 -B4 -O6 -E1 matches BWA-MEM (Sentieon) for concordance.
// ============================================================================
process ALIGN_DEDUP {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)
    tuple path(idx_mbw), path(idx_l2b)
    val(platform)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("${sample_name}.dupblaster.tsv"), emit: stats

    script:
    def prefix = idx_mbw.baseName
    // Split CPUs between the aligner (hot path) and the sorter; dupblaster's own
    // IO threads sit between them and are cheap. Give the sorter ~1/4 of the cores.
    def sort_threads = Math.max(1, (task.cpus as int).intdiv(4))
    def aln_threads  = Math.max(1, (task.cpus as int) - sort_threads)
    // ~40% of the task memory for sort buffers (split across sort threads); the
    // aligner + dupblaster run concurrently in the same pipe and need the rest.
    def sort_mem_mb  = Math.max(768, ((task.memory.toMega() as long) * 0.4 / sort_threads) as int)
    """
    set -euo pipefail
    export TMPDIR=\$PWD
    minibwa-rs map -a -x sr -A1 -B4 -O6 -E1 -Y -t ${aln_threads} \\
        -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:${platform}" \\
        ${prefix} ${reads[0]} ${reads[1]} \\
      | dupblaster --single-end-strategy picard-exact --tmp-dir \$PWD \\
            --stats ${sample_name}.dupblaster.tsv --sample ${sample_name} -o - \\
      | samtools sort -@ ${sort_threads} -m ${sort_mem_mb}M -T ${sample_name}.sort -o ${sample_name}.bam -
    samtools index -@ ${task.cpus} ${sample_name}.bam
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
// PROCESS: RUSTQC_METRICS — single-pass Rust metrics from the deduped BAM
// ============================================================================
process RUSTQC_METRICS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path base_metrics_intervals
    tuple path(reference), path(reference_fai)
    path gc_windows

    output:
    tuple val(sample_name), path("${sample_name}.rustqc.json"), emit: metrics

    script:
    def territory = params.genome_territory ? "--genome-territory ${params.genome_territory}" : ""
    def ref_arg = reference.name != 'NO_FILE' ? "--reference ${reference}" : ""
    def gcwin_arg = gc_windows.name != 'NO_FILE.gcwin' ? "--gc-windows ${gc_windows}" : ""
    """
    bskryb-qc \\
        --bam ${bam} --intervals ${base_metrics_intervals} --sample ${sample_name} \\
        --threads ${task.cpus} \\
        ${territory} ${ref_arg} ${gcwin_arg} \\
        --out ${sample_name}.rustqc.json
    """
}

// ============================================================================
// PROCESS: GINKGO_BINS_TO_PARQUET — SegCopy bin-level CNV -> long-format Parquet
//   (copied from basej-methylqc-rs; "pipeline" field set to basej-dnaqc)
// ============================================================================
process GINKGO_BINS_TO_PARQUET {
    tag "cnv_bins_to_parquet"

    input:
    path(segcopy_file)
    path(raw_counts_merged_file)
    val(bin_size)
    val(read_length)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(user)

    output:
    path("cnv_summary/workspace=${workspace}/workflow_id=${workflow_id}/biosample=*/output.parquet"), emit: cnv_bins_parquet

    script:
    """
    python3 - << 'PYEOF'
import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

workspace = "${workspace}"
workflow_id = "${workflow_id}"
dataset_id = "${dataset_id}"
pipeline_version = "${pipeline_version}"
user = "${user}"
bin_size = int("${bin_size}")
read_length = int("${read_length}")

df_seg = pd.read_csv("${segcopy_file}", sep="\\t")
df_raw = pd.read_csv("${raw_counts_merged_file}", sep="\\t")

sample_cols = [c for c in df_seg.columns if c not in ["CHR", "START", "END"]]
raw_sample_cols = list(df_raw.columns)

print(f"SegCopy has {len(df_seg)} bins and {len(sample_cols)} samples: {sample_cols}")

if len(df_seg) != len(df_raw):
    raise ValueError(f"Bin count mismatch: SegCopy has {len(df_seg)} bins, raw counts has {len(df_raw)} bins")
if len(sample_cols) != len(raw_sample_cols):
    raise ValueError(f"Sample count mismatch: SegCopy has {len(sample_cols)} samples, raw counts has {len(raw_sample_cols)} columns")

for i, sample_col in enumerate(sample_cols):
    biosample = sample_col.replace("_sorted", "")
    raw_col = raw_sample_cols[i]

    df_long = pd.DataFrame({
        "biosample": biosample,
        "dataset_id": dataset_id,
        "pipeline": "basej-dnaqc",
        "pipeline_version": pipeline_version,
        "user": user,
        "bin_size": bin_size,
        "read_length": read_length,
        "chr": df_seg["CHR"],
        "start": df_seg["START"],
        "end": df_seg["END"],
        "ploidy": df_seg[sample_col],
        "raw_count": df_raw[raw_col]
    })

    _str_cols    = ['biosample','dataset_id','pipeline','pipeline_version','user','chr']
    _double_cols = ['ploidy']
    _bigint_cols = ['bin_size','read_length','start','end','raw_count']
    for _col in _double_cols:
        if _col in df_long.columns:
            df_long[_col] = pd.to_numeric(df_long[_col], errors='coerce').astype('float64')
    for _col in _bigint_cols:
        if _col in df_long.columns:
            df_long[_col] = pd.to_numeric(df_long[_col], errors='coerce').astype('Int64')
    for _col in _str_cols:
        if _col in df_long.columns:
            df_long[_col] = df_long[_col].astype('string')
    out_dir = f"cnv_summary/workspace={workspace}/workflow_id={workflow_id}/biosample={biosample}"
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df_long, preserve_index=False),
                   os.path.join(out_dir, "output.parquet"))

print(f"Created CNV summary parquet files for {len(sample_cols)} samples")
PYEOF
    """
}

// ============================================================================
// PROCESS: GINKGO_CNV_SUMMARY — RDS-derived MAPD/SKEW metrics for QC scoring
// ============================================================================
process GINKGO_CNV_SUMMARY {
    tag "ginkgo_cnv_summary"

    input:
    path(ginkgo_rds)

    output:
    path("AllSample-GinkgoSegmentSummary.txt"), emit: summary

    script:
    """
    Rscript /usr/local/bin/cnvSummarizer.R --rds_file ${ginkgo_rds} --out_file AllSample-GinkgoSegmentSummary.txt
    """
}

// ============================================================================
// PROCESS: METRICS_TO_PARQUET — aggregate dnaqc_summary Parquet writer
// ----------------------------------------------------------------------------
// Collect-based rewrite (mirrors basej-methylqc-rs's METHYLQC_METRICS_TO_PARQUET):
// discovers samples from the staged *.rustqc.json and merges the bskryb-qc rustqc
// JSON (align/insert/coverage/chrM/quality-yield/gcbias/preseq/lorenz), dupblaster
// stats (dedup_*), fastp JSON (fastp_*), read counts, and Ginkgo metrics/SegCopy/
// CNV summary (ginkgo_*) into the SAME dnaqc_summary schema (column names +
// types), so the Athena/Iceberg table is unchanged. Also emits the MultiQC
// custom-content table (*_selected_metrics_mqc.txt) and a flat dnaqc_all_metrics.tsv.
// ============================================================================
process METRICS_TO_PARQUET {
    tag "metrics_to_parquet"

    input:
    path(rustqc_jsons)
    path(fastp_jsons)
    path(dupblaster_stats)
    path(read_counts_files)
    path(preseq_files)
    path(ginkgo_metrics)
    path(seg_copy_file, stageAs: 'ginkgo_segcopy')
    path(cnv_summary_file, stageAs: 'ginkgo_cnv_summary')
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(user)

    output:
    path("dnaqc_summary/workspace=*/workflow_id=*/biosample=*/output.parquet"), emit: parquet
    path("*_selected_metrics_mqc.txt"), emit: mqc_metrics
    path("dnaqc_all_metrics.tsv"), emit: summary_tsv
    path("per_biosample_status.csv"), emit: per_biosample_status

    script:
    def has_ginkgo = ginkgo_metrics.name != 'NO_GINKGO_METRICS'
    """
    metrics_to_parquet.py \\
        --dataset-id ${dataset_id} --workspace ${workspace} --workflow-id ${workflow_id} \\
        --pipeline-version ${pipeline_version} --user ${user} \\
        ${has_ginkgo ? "--ginkgo-metrics ${ginkgo_metrics} --segcopy ginkgo_segcopy --cnv-summary ginkgo_cnv_summary" : ""}
    """
}

// ============================================================================
// PROCESS: QC_PLOTS — DNA QC composition + CNV quadrant plots (shared R scripts)
// ----------------------------------------------------------------------------
// Same R scripts and outputs as basej-dnaqc's QC_PLOTS (custom_r_qcplots image):
//   dna_qc_plot.R            -> QC_composition.pdf/_mqc.jpg, DNA-QC_ConsensusScores*,
//                               ConsensusScores_SummaryTableByGroup.txt
//   function_cnv_quadrants_qc.R -> CNV-Quadrants.pdf/_mqc.jpg
// Driven from dnaqc_all_metrics.tsv via prepare_dna_qc_plot_inputs.py (column rename
// only). The Python verdict from METRICS_TO_PARQUET stays authoritative (Parquet,
// per_biosample_status); the R CompositeScore is checked against it in
// qc_score_concordance.tsv. Runs only when Ginkgo ran (the plots need SegCopy + CNV).
// ============================================================================
process QC_PLOTS {
    tag "qc_plots"

    input:
    path(summary_tsv)
    path(seg_copy_file, stageAs: 'ginkgo_segcopy.tsv')
    path(input_csv, stageAs: 'input_samplesheet.csv')
    path(plot_qc_config)

    output:
    path("nf-preseq-pipeline_all_metrics_mqc.txt"), emit: allmetrics_with_cnv
    path("DNA-QC_ConsensusScores.txt"), emit: consensus_scores
    path("DNA-QC_ConsensusScores_SummaryTable_mqc.txt"), emit: consensus_summary
    path("ConsensusScores_SummaryTableByGroup.txt"), emit: consensus_group_summary
    path("QC_composition.pdf"), emit: composition_pdf
    path("QC_composition_mqc.jpg"), emit: composition_jpg
    path("CNV-Quadrants.pdf"), emit: cnv_quadrants_pdf
    path("CNV-Quadrants_mqc.jpg"), emit: cnv_quadrants_jpg
    path("qc_score_concordance.tsv"), emit: concordance

    script:
    """
    set -euo pipefail

    prepare_dna_qc_plot_inputs.py \\
        --summary ${summary_tsv} \\
        --input-csv input_samplesheet.csv \\
        --out-metrics nf-preseq-pipeline_all_metrics_mqc.txt \\
        --out-metadata qc_plot_metadata.csv

    # MAPD/SKEW are already in the metrics table (from the cnvSummarizer output merged by
    # METRICS_TO_PARQUET), so --cnv_summary_file is deliberately NOT passed: with it the
    # R script re-merges the raw summary and discards the SKEW_CNV sentinel fill.
    Rscript /usr/local/bin/dna_qc_plot.R \\
        --seg_copy_file ginkgo_segcopy.tsv \\
        --metrics_file nf-preseq-pipeline_all_metrics_mqc.txt \\
        --metadata_file qc_plot_metadata.csv \\
        --plot_qc_config ${plot_qc_config}

    Rscript /usr/local/bin/function_cnv_quadrants_qc.R \\
        --metrics_file nf-preseq-pipeline_all_metrics_mqc.txt

    compare_qc_scores.py \\
        --summary ${summary_tsv} \\
        --scores DNA-QC_ConsensusScores.txt \\
        --out qc_score_concordance.tsv

    # MultiQC custom-content header so the score distribution renders as a table.
    tmp=\$(mktemp)
    printf '# id: "qc_score_distribution"\\n# section_name: "Total Usable Cells"\\n# description: "Distribution of samples across composite QC score categories. Scores range from 0 (poor quality) to 5 (high quality)."\\n# plot_type: "table"\\n# pconfig:\\n#   id: "qc_score_dist_table"\\n#   title: "Total Usable Cells"\\n' \\
        | cat - DNA-QC_ConsensusScores_SummaryTable_mqc.txt > "\$tmp"
    mv "\$tmp" DNA-QC_ConsensusScores_SummaryTable_mqc.txt
    """
}

// ============================================================================
// PROCESS: PRESEQ — library complexity via the actual preseq tool
//   bskryb-qc's gc_extrap port is faithful GIVEN the same reads, but the minibwa
//   BAM (-a secondary/supplementary + dupblaster-marked, not removed, dups) is a
//   different read set than the Sentieon rmdup BAM the reference dnaqc feeds preseq,
//   so preseq_count undershoots. Run the real preseq tool on a cleaned primary BAM
//   for a like-for-like library-complexity number.
// ============================================================================
process PRESEQ {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    tuple val(sample_name), path("${sample_name}_preseq.txt"), emit: complexity

    script:
    """
    set +u

    echo "0" > ${sample_name}_preseq.txt

    # Keep primary, mapped, PROPERLY-PAIRED, non-duplicate reads (-f 2 -F 3844) so the
    # input matches the Sentieon rmdup primary BAM. -f 2 drops discordant/large-insert
    # pairs that bam2mr merges into fragments wider than gc_extrap's max_width
    # ("read of width N max_width set too small" crash); then sort the .mr. -w is
    # belt-and-suspenders for any residual wide fragment (seg_len is 100kb).
    # Single-end input (Ultima CRAM path) has no proper-pair flag, so -f 2 would drop
    # every read; there the filter is primary, mapped, non-duplicate only (-F 3844).
    if [ "\$(samtools view -c -f 1 -F 0x900 ${bam})" -gt 0 ]; then
        samtools view -b -f 2 -F 3844 ${bam} > ${sample_name}.primary.bam
    else
        samtools view -b -F 3844 ${bam} > ${sample_name}.primary.bam
    fi

    {
        bam2mr -seg_len 100000 -o ${sample_name}.mr.unsorted ${sample_name}.primary.bam 2>/dev/null && \\
        LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k6,6 ${sample_name}.mr.unsorted > ${sample_name}.mr && \\
        preseq gc_extrap -w 300000 -o ${sample_name}.curve ${sample_name}.mr && \\
        tail -n 1 ${sample_name}.curve | cut -f 2 | awk '{printf "%.0f\\n", \$1}' > ${sample_name}_preseq.txt
    } || {
        echo "Preseq failed (likely low complexity library), using 0"
    }

    echo "Estimated library complexity: \$(cat ${sample_name}_preseq.txt)"
    """
}

// ============================================================================
// PROCESS: MULTIQC — QC report (title: basej-dnaqc)
//   (copied from basej-methylqc-rs; title/section adjusted to basej-dnaqc)
// ============================================================================
process MULTIQC {
    tag "multiqc"

    input:
    path(mqc_files)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(nextflow_version)
    val(architecture)
    val(run_ginkgo)
    val(seqkit_container)
    val(fastp_container)
    val(align_container)
    val(rustqc_container)
    val(metrics_container)
    path(logo)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data
    path("tool_mqc_versions.yml"), emit: versions

    script:
    """
    set -euo pipefail

    MULTIQC_VERSION=\$(multiqc --version 2>&1 | head -n 1 | sed -E 's/^multiqc, version[[:space:]]+//' | tr -d '\r')
    test -n "\$MULTIQC_VERSION"

    # Legacy BioSkryb/MultiQC convention: grouped tool versions are auto-discovered
    # and rendered in the standard Software Versions section.
    cat > tool_mqc_versions.yml << EOF
basej-dnaqc:
  pipeline: "${pipeline_version}"
  nextflow: "${nextflow_version}"
  seqkit: "2.13.0"
  fastp: "0.20.1"
  minibwa-rs: "0.6-r416"
  dupblaster: "0.1.1"
  samtools-align: "1.16.1"
  bskryb-qc-package: "0.1.5"
  bskryb-rustqc-image-release: "0.1.5"
  preseq: "2.0.3"
  bam2mr: "preseq-2.0.3"
  samtools-preseq: "1.9"
  multiqc: "\$MULTIQC_VERSION"
  execution-architecture: "${architecture}"
  seqkit-container: "${seqkit_container}"
  fastp-container: "${fastp_container}"
  align-container: "${align_container}"
  rustqc-container: "${rustqc_container}"
  preseq-container: "basejumper_preseq_bam2mr_0.1"
  metrics-container: "${metrics_container}"
  multiqc-container: "quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"
EOF

    if [ "${run_ginkgo}" = "true" ]; then
        cat >> tool_mqc_versions.yml << EOF
  bedtools: "2.28.0"
  ginkgo: "0.0.2"
  ginkgo-core-container: "basejumper_ginkgo_0.3.1"
  ginkgo-parser-container: "basejumper_ginko_parser_0.2.1"
  ginkgo-rds-metrics-container: "basejumper_ginkgo_0.3.1"
  ginkgo-summary-container: "basejumper_custom_r_qcplots_0.3.0"
EOF
    fi

    cat > multiqc_config.yaml << EOF
custom_logo_title: 'BioSkryb Genomics'
custom_logo: ${logo}
custom_logo_width: 260

title: "basej-dnaqc v${pipeline_version}"
report_header_info:
  - Dataset ID: "${dataset_id}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false
skip_generalstats: true

fn_clean_exts:
  - ".fastq.gz"
  - ".fq.gz"
  - ".bam"
  - "_fastp"
  - "_trim"
  - ".txt"
  - ".json"

module_order:
  - custom_content
  - fastp

report_section_order:
  dnaqc_summary:
    order: 1000
  QC_composition_mqc.jpg:
    order: 900
  qc_score_distribution:
    order: 850
  CNV-Quadrants_mqc.jpg:
    order: 800

custom_data:
  QC_composition_mqc.jpg:
    section_name: "QC Composition"
    description: "Per-sample QC metric distributions grouped by QC tier, with the Ginkgo copy-number heatmap."
  CNV-Quadrants_mqc.jpg:
    section_name: "CNV Quadrants"
    description: "MAPD vs. segment unevenness (skew) of Ginkgo CNV bins per cell; shaded regions mark the QC gates."

table_cond_formatting_rules:
  QC_Status:
    pass:
      - s_eq: "PASS"
    warn:
      - s_eq: "Borderline"
    fail:
      - s_eq: "FAIL"
EOF

    multiqc . -n multiqc_report.html -c multiqc_config.yaml --force
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================
workflow {
    main:
    def resolved_pipeline_version = params.pipeline_version ?: workflow.manifest.version
    def resolved_nextflow_version = workflow.nextflow.version.toString()
    def seqkit_container = "basejumper_seqkit-2.13.0"
    def fastp_container = "quay.io/biocontainers/fastp:0.20.1--h8b12597_0"
    def align_container = "basejumper_minibwa-rustqc_0.1.5"
    def rustqc_container = "basejumper_bskryb-rustqc-x86_0.1.5"
    def metrics_container = "basejumper_custom_parabricks-metrics_1.0.3"

    if (!params.input_csv) {
        error "ERROR: --input_csv is required"
    }

    // Platform is auto-detected from the CSV (same rule as basej-dnaqc): any row with a
    // non-empty `cram` column makes this an Ultima run (biosampleName,cram[,crai]);
    // otherwise Illumina/Element FASTQ (biosampleName,read1,read2). One platform per run.
    def csv_rows = file(params.input_csv, checkIfExists: true).splitCsv(header: true)
    def is_ultima = csv_rows.any { row -> row.cram?.trim() }
    log.info "basej-dnaqc: auto-detected platform: ${is_ultima ? 'Ultima (CRAM)' : 'Illumina (FASTQ)'}"

    ch_intervals = Channel.value(file(params.intervals))

    if (is_ultima) {
        // ===== Ultima: pre-aligned CRAM -> subsample -> decode + re-dedup (no realign) =====
        ch_cram = Channel
            .fromPath(params.input_csv, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row.cram?.trim() }
            .map { row ->
                def name = (row.biosampleName ?: row.biosample)?.trim()
                if (!name) {
                    error "ERROR: input_csv row is missing biosampleName: ${row}"
                }
                def cram = file(row.cram.trim())
                def crai = row.crai?.trim() ? file(row.crai.trim()) : file("${row.cram.trim()}.crai")
                tuple(name, cram, crai)
            }
        ch_cram_ref = Channel.value([
            file(params.reference_fasta, checkIfExists: true),
            file("${params.reference_fasta}.fai", checkIfExists: true)
        ])

        SAMTOOLS_SUBSAMPLE_CRAM(
            ch_cram.map { s, c, i -> tuple(s, c, i, params.n_reads) },
            params.samtools_seed,
            ch_cram_ref
        )
        CRAM_DEDUP(SAMTOOLS_SUBSAMPLE_CRAM.out.cram, ch_cram_ref)

        ch_bam          = CRAM_DEDUP.out.bam
        ch_dedup_stats  = CRAM_DEDUP.out.stats
        ch_read_counts  = SAMTOOLS_SUBSAMPLE_CRAM.out.read_counts
        ch_fastp_json   = channel.empty()   // no FASTQ -> no fastp; Q30 comes from QualityYield
    } else {
        // ===== Illumina/Element: FASTQ -> subsample -> fastp -> align|dedup|sort =====
        // Input CSV columns: biosampleName (or biosample), read1, read2
        //   Multi-lane: read1/read2 may hold pipe-delimited ("|") lane paths, e.g.
        //     S1,s3://b/L001_R1.fq.gz|s3://b/L002_R1.fq.gz,s3://b/L001_R2.fq.gz|s3://b/L002_R2.fq.gz
        ch_raw = Channel
            .fromPath(params.input_csv, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row.read1?.trim() && row.read2?.trim() }
            .map { row ->
                def name = (row.biosampleName ?: row.biosample)?.trim()
                if (!name) {
                    error "ERROR: input_csv row is missing biosampleName: ${row}"
                }
                def r1 = row.read1.tokenize('|').collect { file(it.trim()) }
                def r2 = row.read2.tokenize('|').collect { file(it.trim()) }
                if (r1.size() != r2.size()) {
                    error "ERROR: biosample '${name}' has ${r1.size()} read1 lane(s) but ${r2.size()} read2 lane(s)"
                }
                tuple(name, r1, r2)
            }

        ch_raw_branched = ch_raw.branch { name, r1, r2 ->
            multilane: r1.size() > 1
            singlelane: true
        }
        MERGE_MULTILANE_FASTQ(ch_raw_branched.multilane)
        ch_reads = ch_raw_branched.singlelane
            .map { name, r1, r2 -> tuple(name, [r1[0], r2[0]]) }
            .mix(MERGE_MULTILANE_FASTQ.out.reads)

        ch_minibwa_index = Channel.value([file("${params.minibwa_index}.mbw"), file("${params.minibwa_index}.l2b")])

        // Subsample -> trim -> (align | dedup | sort in one streaming pass)
        SEQKIT_SAMPLE(ch_reads.map { s, r -> tuple(s, r, params.n_reads) }, params.seqkit_sample_seed)
        FASTP_TRIM(SEQKIT_SAMPLE.out.reads)
        ALIGN_DEDUP(FASTP_TRIM.out.reads, ch_minibwa_index, params.platform ?: 'ILLUMINA')

        ch_bam          = ALIGN_DEDUP.out.bam
        ch_dedup_stats  = ALIGN_DEDUP.out.stats
        ch_read_counts  = SEQKIT_SAMPLE.out.read_counts
        ch_fastp_json   = FASTP_TRIM.out.json
    }

    // GcBias reference + windows index (opt-in)
    ch_reference = params.run_gcbias
        ? Channel.value([file(params.reference_fasta), file("${params.reference_fasta}.fai")])
        : Channel.value([file("${projectDir}/assets/NO_FILE"), file("${projectDir}/assets/NO_FILE.fai")])
    if (params.run_gcbias) {
        if (params.gc_windows_index) {
            ch_gcwin = Channel.value(file(params.gc_windows_index))
        } else {
            GCBIAS_INDEX(ch_reference, ch_intervals)
            ch_gcwin = GCBIAS_INDEX.out.index.first()
        }
    } else {
        ch_gcwin = Channel.value(file("${projectDir}/assets/NO_FILE.gcwin"))
    }

    RUSTQC_METRICS(ch_bam, ch_intervals, ch_reference, ch_gcwin)

    // Library complexity via the actual preseq tool (not the bskryb-qc port)
    PRESEQ(ch_bam)

    // ---- Ginkgo CNV (runs for GRCh38, GRCm39, ARSUCD2 when !skip_cnv) ----
    ch_ginkgo_metrics_qc = Channel.value(file("${projectDir}/assets/NO_GINKGO_METRICS"))
    ch_ginkgo_segcopy_qc = Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_ginkgo_cnv_summary_qc = Channel.value(file("${projectDir}/assets/NO_FILE"))

    // Publishable CNV channels — empty unless Ginkgo runs for this genome.
    ch_ginkgo_rds = channel.empty()
    ch_ginkgo_segcopy = channel.empty()
    ch_ginkgo_cnv_plots = channel.empty()
    ch_cnv_summary_parquet = channel.empty()

    def effective_skip_cnv = params.skip_cnv != null ? params.skip_cnv : (params.genome == "pUC19")
    def run_ginkgo = !effective_skip_cnv && params.genome in ["GRCh38", "GRCm39", "ARSUCD2"]

    if (run_ginkgo) {
        ch_binref    = channel.fromPath(params.ginko_ref_dir + "variable_" + params.bin_size + "_" + params.read_length + "_bwa")
        ch_gcref     = channel.fromPath(params.ginko_ref_dir + "GC_variable_" + params.bin_size + "_" + params.read_length + "_bwa")
        ch_boundsref = channel.fromPath(params.ginko_ref_dir + "bounds_variable_" + params.bin_size + "_" + params.read_length + "_bwa")

        GINKO_NOPUBLISH(
            ch_bam,
            params.bin_size,
            ch_binref,
            ch_gcref,
            ch_boundsref,
            params.min_ploidy,
            params.max_ploidy,
            params.min_bin_width,
            params.is_haplotype
        )

        ch_ginkgo_metrics_qc = GINKO_NOPUBLISH.out.metrics
        ch_ginkgo_segcopy_qc = GINKO_NOPUBLISH.out.segcopy

        ch_ginkgo_rds = GINKO_NOPUBLISH.out.rds
        ch_ginkgo_segcopy = GINKO_NOPUBLISH.out.segcopy
        ch_ginkgo_cnv_plots = GINKO_NOPUBLISH.out.graph

        GINKGO_CNV_SUMMARY(GINKO_NOPUBLISH.out.rds)
        ch_ginkgo_cnv_summary_qc = GINKGO_CNV_SUMMARY.out.summary

        GINKGO_BINS_TO_PARQUET(
            GINKO_NOPUBLISH.out.segcopy,
            GINKO_NOPUBLISH.out.raw_counts_merged,
            params.bin_size,
            params.read_length,
            params.dataset_id,
            params.workspace,
            params.workflow_id,
            resolved_pipeline_version,
            params.pipeline_user
        )

        ch_cnv_summary_parquet = GINKGO_BINS_TO_PARQUET.out.cnv_bins_parquet
    }

    // ---- dnaqc_summary Parquet (aggregate, collect-based) ----
    METRICS_TO_PARQUET(
        RUSTQC_METRICS.out.metrics.map { s, j -> j }.collect(),
        ch_fastp_json.map { s, j -> j }.collect().ifEmpty([]),
        ch_dedup_stats.map { s, t -> t }.collect(),
        ch_read_counts.map { s, r -> r }.collect(),
        PRESEQ.out.complexity.map { s, t -> t }.collect(),
        ch_ginkgo_metrics_qc,
        ch_ginkgo_segcopy_qc,
        ch_ginkgo_cnv_summary_qc,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        resolved_pipeline_version,
        params.pipeline_user
    )

    // ---- QC plots (composition / consensus scores / CNV quadrants) ----
    // Needs Ginkgo SegCopy + CNV summary, so it runs only when Ginkgo ran.
    ch_qc_plots = channel.empty()
    ch_qc_plots_mqc = channel.empty()
    if (run_ginkgo) {
        QC_PLOTS(
            METRICS_TO_PARQUET.out.summary_tsv,
            GINKO_NOPUBLISH.out.segcopy,
            file(params.input_csv, checkIfExists: true),
            file("${projectDir}/assets/plot_qc_config.json", checkIfExists: true)
        )
        ch_qc_plots = QC_PLOTS.out.allmetrics_with_cnv
            .mix(QC_PLOTS.out.consensus_scores)
            .mix(QC_PLOTS.out.consensus_summary)
            .mix(QC_PLOTS.out.consensus_group_summary)
            .mix(QC_PLOTS.out.composition_pdf)
            .mix(QC_PLOTS.out.composition_jpg)
            .mix(QC_PLOTS.out.cnv_quadrants_pdf)
            .mix(QC_PLOTS.out.cnv_quadrants_jpg)
            .mix(QC_PLOTS.out.concordance)
        ch_qc_plots_mqc = QC_PLOTS.out.composition_jpg
            .mix(QC_PLOTS.out.cnv_quadrants_jpg)
            .mix(QC_PLOTS.out.consensus_summary)
    }

    // ---- MultiQC ----
    ch_mqc_inputs = ch_fastp_json.map { s, j -> j }.collect().ifEmpty([])
        .mix(METRICS_TO_PARQUET.out.mqc_metrics.collect())
        .mix(ch_qc_plots_mqc)
        .collect()

    MULTIQC(
        ch_mqc_inputs,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        resolved_pipeline_version,
        resolved_nextflow_version,
        params.architecture,
        run_ginkgo,
        seqkit_container,
        fastp_container,
        align_container,
        rustqc_container,
        metrics_container,
        file("${projectDir}/assets/bioskryb_logo-tagline.png", checkIfExists: true)
    )

    publish:
    bam_files = ch_bam
        .map { sample_name, bam, bai -> [biosampleName: sample_name, bam: bam, bai: bai] }
    rustqc_metrics = RUSTQC_METRICS.out.metrics
        .map { sample_name, metrics -> [biosampleName: sample_name, metrics: metrics] }
    dedup_metrics = ch_dedup_stats
        .map { sample_name, metrics -> [biosampleName: sample_name, metrics: metrics] }
    preseq_complexity = PRESEQ.out.complexity
        .map { sample_name, preseq -> [biosampleName: sample_name, preseq: preseq] }
    dnaqc_summary = METRICS_TO_PARQUET.out.parquet
    dnaqc_summary_tsv = METRICS_TO_PARQUET.out.summary_tsv
    cnv_summary = ch_cnv_summary_parquet
    ginkgo_rds = ch_ginkgo_rds
    ginkgo_segcopy = ch_ginkgo_segcopy
    cnv_plots_per_cell = ch_ginkgo_cnv_plots
    per_biosample_status = METRICS_TO_PARQUET.out.per_biosample_status
    qc_plots = ch_qc_plots
    software_versions = MULTIQC.out.versions
    multiqc_report = MULTIQC.out.report
}

// ============================================================================
// OUTPUT CONFIGURATION
// ============================================================================
// Layout mirrors basej-dnaqc exactly so the platform (BaseJumper connector +
// lakehouse) reads this pipeline with no special-casing:
//   tables/<table>/workspace=/workflow_id=/biosample=/  -> Glue/Iceberg ingestion
//   workflow_outputs/<ws>/<wf>/index/bam.csv            -> downstream BAM wiring
//   workflow_outputs/<ws>/<wf>/reports/                  -> MultiQC report
//   workflow_outputs/<ws>/<wf>/cnv/ginkgo_cnv_plots/     -> declared CNV artifact
// ============================================================================
output {
    bam_files {
        path "bam/${params.workspace}/dna/tool=minibwa-rs/pipeline=dnaqc"
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

    preseq_complexity {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/preseq"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "preseq",
             tool: "preseq"
    }

    // Per-biosample dnaqc_summary parquet -> Athena/Iceberg via the Glue job.
    // MUST land at the bucket-root `tables/` prefix (same as basej-dnaqc).
    dnaqc_summary {
        path "tables"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "dnaqc_summary",
             reference: params.genome
    }

    // Flat TSV mirroring the parquet schema (nf-test / parity validation)
    dnaqc_summary_tsv {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/qc_metrics"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "dnaqc_summary_tsv"
    }

    cnv_summary {
        path "tables"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_summary",
             tool: "ginkgo",
             reference: params.genome
    }

    ginkgo_rds {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/cnv/ginkgo_rds"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_rds",
             tool: "ginkgo"
    }

    ginkgo_segcopy {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/cnv/ginkgo_segcopy"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_segcopy",
             tool: "ginkgo"
    }

    // Per-cell CNV profile plots (tar.gz of one JPEG per biosample)
    cnv_plots_per_cell {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/cnv/ginkgo_cnv_plots"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_plots",
             tool: "ginkgo"
    }

    // Composition / consensus-score / CNV-quadrant plots (same layout as basej-dnaqc)
    qc_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/qc_plots"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "qc_plots"
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
