// ============================================================================
// BASEJ-LINEAGE LOCAL MODULES  (matrix-consumer pipeline)
// ============================================================================
// All processes are defined inline (no publishDir — outputs are handled via
// Seqera Platform publish blocks in main.nf).
// Containers and resources are configured in nextflow.config under withName
// directives. The publish_dir / enable_publish / params.timestamp machinery
// from the source subworkflows is intentionally stripped per basej convention.
//
// Two branches merged from Isai's standalone subworkflows:
//   * lineage_from_matrices     — phylogeny / lineage from NR/NV matrices
//   * mutsignatures_from_matrices — COSMIC signatures from a binary matrix
// ============================================================================

nextflow.enable.dsl=2

// ============================================================================
// LINEAGE BRANCH
// ============================================================================

// ---------------------------------------------------------------------------
// GENOTYPE_TABLE_FROM_ANNOTATED_VCF
// Emit a two-column TSV (VARIANT_ID<tab>GT) per sample from an annotated VCF.
// VARIANT_ID = CHROM_POS_REF_ALT. One task per sample. Consumed by POSTPROCESS.
// ---------------------------------------------------------------------------
process GENOTYPE_TABLE_FROM_ANNOTATED_VCF {
    tag "${group}:${sample_name}"

    input:
    tuple val(group), val(sample_name), path(vcf), path(tbi)

    output:
    tuple val(group), val(sample_name),
          path("${sample_name}_genotype.tsv"),
          emit: genotype_table

    script:
    """
    set -euo pipefail

    printf 'VARIANT_ID\tGT\n' > ${sample_name}_genotype.tsv
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' ${vcf} \
        | awk 'BEGIN{OFS="\t"} {print \$1"_"\$2"_"\$3"_"\$4, \$5}' \
        >> ${sample_name}_genotype.tsv
    """
}

// ---------------------------------------------------------------------------
// SEQUOIA_PHYLOGENY_SNV
// Build phylogeny for one NR/NV matrix pair (single label) after subsetting to
// SNV rows. genotype_bin (5th tuple element) optionally supplies pre-built calls.
// ---------------------------------------------------------------------------
process SEQUOIA_PHYLOGENY_SNV {
    tag "${group}_${label}"

    input:
    tuple val(group), val(label), path(nr_matrix), path(nv_matrix), path(genotype_bin_matrix)
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(create_multi_tree)
    val(mpboot_path)
    val(mpboot_bootstrap)

    output:
    tuple val(group), path("output_snv_${label}"), emit: phylogeny_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Always create the output dir so the Nextflow output declaration always resolves
    mkdir -p output_snv_${label}
    touch output_snv_${label}/no_results

    # ── Subset to SNV rows (len(REF)==1 && len(ALT)==1) ──────────────────────
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)==1 && length(alt)==1) print }' "${nr_matrix}" > snv_nr.tsv
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)==1 && length(alt)==1) print }' "${nv_matrix}" > snv_nv.tsv

    n_snv=\$(tail -n +2 snv_nr.tsv | wc -l)
    echo "[snv:${label}] SNV variants: \${n_snv}"

    if [ "\${n_snv}" -eq 0 ]; then
        echo "[snv:${label}] No SNVs — skipping"
        exit 0
    fi

    # ── Optional pre-built genotype_bin (non-empty sentinel file) ────────────
    genotype_bin_arg=""
    if [ -s "${genotype_bin_matrix}" ]; then
        awk 'NR==1 { print; next }
             { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
               if (length(ref)==1 && length(alt)==1) print }' "${genotype_bin_matrix}" > snv_genotype_bin.tsv
        genotype_bin_arg="--input_genotype_bin snv_genotype_bin.tsv"
    fi

    Rscript /usr/local/bin/rscript_sequoia_build_phylogeny_only.R \\
        --donor_id        "${group}_${label}" \\
        --input_nr        snv_nr.tsv \\
        --input_nv        snv_nv.tsv \\
        --output_dir      output_snv_${label}/ \\
        --only_snvs       TRUE \\
        --gender          "${gender}" \\
        --vaf_absent      ${vaf_absent} \\
        --vaf_present     ${vaf_present} \\
        --create_multi_tree           ${create_multi_tree} \\
        --mpboot_path     "${mpboot_path}" \\
        --mpboot_bootstrap            ${mpboot_bootstrap} \\
        \${genotype_bin_arg}
    """
}

// ---------------------------------------------------------------------------
// SEQUOIA_PHYLOGENY_INDEL
// As SEQUOIA_PHYLOGENY_SNV but subset to indel rows (len(REF)>1 || len(ALT)>1).
// ---------------------------------------------------------------------------
process SEQUOIA_PHYLOGENY_INDEL {
    tag "${group}_${label}"

    input:
    tuple val(group), val(label), path(nr_matrix), path(nv_matrix), path(genotype_bin_matrix)
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(create_multi_tree)
    val(mpboot_path)
    val(mpboot_bootstrap)

    output:
    tuple val(group), path("output_indel_${label}"), emit: phylogeny_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Always create the output dir so the Nextflow output declaration always resolves
    mkdir -p output_indel_${label}
    touch output_indel_${label}/no_results

    # ── Subset to indel rows (len(REF)>1 || len(ALT)>1) ─────────────────────
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)>1 || length(alt)>1) print }' "${nr_matrix}" > indel_nr.tsv
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)>1 || length(alt)>1) print }' "${nv_matrix}" > indel_nv.tsv

    n_indel=\$(tail -n +2 indel_nr.tsv | wc -l)
    echo "[indel:${label}] Indel variants: \${n_indel}"

    if [ "\${n_indel}" -eq 0 ]; then
        echo "[indel:${label}] No indels — skipping"
        exit 0
    fi

    # ── Optional pre-built genotype_bin (non-empty sentinel file) ────────────
    genotype_bin_arg=""
    if [ -s "${genotype_bin_matrix}" ]; then
        awk 'NR==1 { print; next }
             { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
               if (length(ref)>1 || length(alt)>1) print }' "${genotype_bin_matrix}" > indel_genotype_bin.tsv
        genotype_bin_arg="--input_genotype_bin indel_genotype_bin.tsv"
    fi

    Rscript /usr/local/bin/rscript_sequoia_build_phylogeny_only.R \\
        --donor_id        "${group}_${label}" \\
        --input_nr        indel_nr.tsv \\
        --input_nv        indel_nv.tsv \\
        --output_dir      output_indel_${label}/ \\
        --only_snvs       FALSE \\
        --gender          "${gender}" \\
        --vaf_absent      ${vaf_absent} \\
        --vaf_present     ${vaf_present} \\
        --create_multi_tree           ${create_multi_tree} \\
        --mpboot_path     "${mpboot_path}" \\
        --mpboot_bootstrap            ${mpboot_bootstrap} \\
        \${genotype_bin_arg}
    """
}

// ---------------------------------------------------------------------------
// SEQUOIA_PHYLOGENY_BOTH
// As SEQUOIA_PHYLOGENY_SNV but uses all variants (SNVs + indels).
// ---------------------------------------------------------------------------
process SEQUOIA_PHYLOGENY_BOTH {
    tag "${group}_${label}"

    input:
    tuple val(group), val(label), path(nr_matrix), path(nv_matrix), path(genotype_bin_matrix)
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(create_multi_tree)
    val(mpboot_path)
    val(mpboot_bootstrap)

    output:
    tuple val(group), path("output_both_${label}"), emit: phylogeny_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Always create the output dir so the Nextflow output declaration always resolves
    mkdir -p output_both_${label}
    touch output_both_${label}/no_results

    n_vars=\$(tail -n +2 "${nr_matrix}" | wc -l)
    echo "[both:${label}] Total variants: \${n_vars}"

    if [ "\${n_vars}" -eq 0 ]; then
        echo "[both:${label}] No variants — skipping"
        exit 0
    fi

    # ── Optional pre-built genotype_bin (non-empty sentinel file) ────────────
    genotype_bin_arg=""
    if [ -s "${genotype_bin_matrix}" ]; then
        genotype_bin_arg="--input_genotype_bin ${genotype_bin_matrix}"
    fi

    Rscript /usr/local/bin/rscript_sequoia_build_phylogeny_only.R \\
        --donor_id        "${group}_${label}" \\
        --input_nr        "${nr_matrix}" \\
        --input_nv        "${nv_matrix}" \\
        --output_dir      output_both_${label}/ \\
        --only_snvs       FALSE \\
        --gender          "${gender}" \\
        --vaf_absent      ${vaf_absent} \\
        --vaf_present     ${vaf_present} \\
        --create_multi_tree           ${create_multi_tree} \\
        --mpboot_path     "${mpboot_path}" \\
        --mpboot_bootstrap            ${mpboot_bootstrap} \\
        \${genotype_bin_arg}
    """
}

// ---------------------------------------------------------------------------
// SEQUOIA_VARIANT_PLACEMENT_SNV
// Place SNV variants on the phylogeny tree via rscript_minimal_variant_placement.R.
// Extra val() inputs are kept for workflow compatibility (ignored by the R script).
// ---------------------------------------------------------------------------
process SEQUOIA_VARIANT_PLACEMENT_SNV {
    tag "${group}"

    input:
    // stageAs: distinct names so sentinel file('/dev/null') for tree + phylo does not collide (both basename "null").
    tuple val(group),
          path(tree_file, stageAs: 'placement_input_tree'),
          path(phylo_placed_variants_all, stageAs: 'placement_phylo_placed_variants_all.tsv'),
          path(nr_matrix, stageAs: 'placement_nr_matrix.tsv'),
          path(nv_matrix, stageAs: 'placement_nv_matrix.tsv')
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(tree_mut_pval)
    val(keep_ancestral)
    val(create_multi_tree)
    val(genotype_conv_prob)
    val(min_pval_for_true_somatic)
    val(min_variant_reads_shared)
    val(min_vaf_shared)

    output:
    tuple val(group), path("output_snv_placement_*"), emit: placement_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob always resolves even when the run is skipped
    touch output_snv_placement_no_results

    # ── Guard: tree file must be non-empty ────────────────────────────────────
    tree_bytes=\$(wc -c < "${tree_file}" || echo 0)
    if [ "\${tree_bytes}" -eq 0 ]; then
        echo "[snv_placement] Tree file is empty (0 bytes) — skipping"
        exit 0
    fi

    NR_IN="${nr_matrix}"
    NV_IN="${nv_matrix}"
    if [ -s "${phylo_placed_variants_all}" ]; then
        tail -n +2 "${phylo_placed_variants_all}" | tr -d '\r' | awk -F'\\t' 'NF>=4 { print \$1"_"\$2"_"\$3"_"\$4 }' | sort -u > exclude_phylo_variant_ids.txt
        awk -F'\\t' 'NR==FNR{ex[\$1]=1;next} NR==1{print;next} !(\$1 in ex)' exclude_phylo_variant_ids.txt "${nr_matrix}" > nr_for_placement.tsv
        awk -F'\\t' 'NR==FNR{ex[\$1]=1;next} NR==1{print;next} !(\$1 in ex)' exclude_phylo_variant_ids.txt "${nv_matrix}" > nv_for_placement.tsv
        NR_IN="nr_for_placement.tsv"
        NV_IN="nv_for_placement.tsv"
        echo "[snv_placement] Excluded \$(wc -l < exclude_phylo_variant_ids.txt) distinct variant IDs present in phylo *_placed_variants_all.tsv"
    fi

    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)==1 && length(alt)==1) print }' "\${NR_IN}" > nr_run.tsv
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)==1 && length(alt)==1) print }' "\${NV_IN}" > nv_run.tsv

    n_snv=\$(tail -n +2 nr_run.tsv | wc -l)
    echo "[snv_placement] SNV variants in matrix: \${n_snv}"

    if [ "\${n_snv}" -eq 0 ]; then
        echo "[snv_placement] No SNVs — skipping"
        exit 0
    fi

    mkdir -p output_snv_placement_${group}

    # Input tree from SEQUOIA_PHYLOGENY_* (branch-length Newick or MPBoot .treefile fallback)
    cp "${tree_file}" "output_snv_placement_${group}/"

    # Phylogeny placement table from the same scheme as tree_file (filtered matrix, treemut on built tree)
    if [ -s "${phylo_placed_variants_all}" ]; then
        cp "${phylo_placed_variants_all}" "output_snv_placement_${group}/"
    fi

    Rscript /usr/local/bin/rscript_minimal_variant_placement.R \\
        --donor_id                   "${group}" \\
        --input_nr                   nr_run.tsv \\
        --input_nv                   nv_run.tsv \\
        --input_tree                 "${tree_file}" \\
        --output_dir                 "output_snv_placement_${group}/" \\
        --variant_type               snv \\
        --tree_mut_pval              ${tree_mut_pval} \\
        --phylo_placed_variants_tsv  "${phylo_placed_variants_all}"
    """
}

// ---------------------------------------------------------------------------
// SEQUOIA_VARIANT_PLACEMENT_INDEL
// ---------------------------------------------------------------------------
process SEQUOIA_VARIANT_PLACEMENT_INDEL {
    tag "${group}"

    input:
    // stageAs: distinct names so sentinel file('/dev/null') for tree + phylo does not collide (both basename "null").
    tuple val(group),
          path(tree_file, stageAs: 'placement_input_tree'),
          path(phylo_placed_variants_all, stageAs: 'placement_phylo_placed_variants_all.tsv'),
          path(nr_matrix, stageAs: 'placement_nr_matrix.tsv'),
          path(nv_matrix, stageAs: 'placement_nv_matrix.tsv')
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(tree_mut_pval)
    val(keep_ancestral)
    val(create_multi_tree)
    val(genotype_conv_prob)
    val(min_pval_for_true_somatic)
    val(min_variant_reads_shared)
    val(min_vaf_shared)

    output:
    tuple val(group), path("output_indel_placement_*"), emit: placement_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob always resolves even when the run is skipped
    touch output_indel_placement_no_results

    # ── Guard: tree file must be non-empty ────────────────────────────────────
    tree_bytes=\$(wc -c < "${tree_file}" || echo 0)
    if [ "\${tree_bytes}" -eq 0 ]; then
        echo "[indel_placement] Tree file is empty (0 bytes) — skipping"
        exit 0
    fi

    NR_IN="${nr_matrix}"
    NV_IN="${nv_matrix}"
    if [ -s "${phylo_placed_variants_all}" ]; then
        tail -n +2 "${phylo_placed_variants_all}" | tr -d '\r' | awk -F'\\t' 'NF>=4 { print \$1"_"\$2"_"\$3"_"\$4 }' | sort -u > exclude_phylo_variant_ids.txt
        awk -F'\\t' 'NR==FNR{ex[\$1]=1;next} NR==1{print;next} !(\$1 in ex)' exclude_phylo_variant_ids.txt "${nr_matrix}" > nr_for_placement.tsv
        awk -F'\\t' 'NR==FNR{ex[\$1]=1;next} NR==1{print;next} !(\$1 in ex)' exclude_phylo_variant_ids.txt "${nv_matrix}" > nv_for_placement.tsv
        NR_IN="nr_for_placement.tsv"
        NV_IN="nv_for_placement.tsv"
        echo "[indel_placement] Excluded \$(wc -l < exclude_phylo_variant_ids.txt) distinct variant IDs present in phylo *_placed_variants_all.tsv"
    fi

    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)>1 || length(alt)>1) print }' "\${NR_IN}" > nr_run.tsv
    awk 'NR==1 { print; next }
         { n=split(\$1,a,"_"); ref=a[n-1]; alt=a[n];
           if (length(ref)>1 || length(alt)>1) print }' "\${NV_IN}" > nv_run.tsv

    n_indel=\$(tail -n +2 nr_run.tsv | wc -l)
    echo "[indel_placement] INDEL variants in matrix: \${n_indel}"

    if [ "\${n_indel}" -eq 0 ]; then
        echo "[indel_placement] No INDELs — skipping"
        exit 0
    fi

    mkdir -p output_indel_placement_${group}

    cp "${tree_file}" "output_indel_placement_${group}/"

    if [ -s "${phylo_placed_variants_all}" ]; then
        cp "${phylo_placed_variants_all}" "output_indel_placement_${group}/"
    fi

    Rscript /usr/local/bin/rscript_minimal_variant_placement.R \\
        --donor_id                   "${group}" \\
        --input_nr                   nr_run.tsv \\
        --input_nv                   nv_run.tsv \\
        --input_tree                 "${tree_file}" \\
        --output_dir                 "output_indel_placement_${group}/" \\
        --variant_type               indel \\
        --tree_mut_pval              ${tree_mut_pval} \\
        --phylo_placed_variants_tsv  "${phylo_placed_variants_all}"
    """
}

// ---------------------------------------------------------------------------
// SEQUOIA_VARIANT_PLACEMENT_BOTH
// ---------------------------------------------------------------------------
process SEQUOIA_VARIANT_PLACEMENT_BOTH {
    tag "${group}"

    input:
    // stageAs: distinct names so sentinel file('/dev/null') for tree + phylo does not collide (both basename "null").
    tuple val(group),
          path(tree_file, stageAs: 'placement_input_tree'),
          path(phylo_placed_variants_all, stageAs: 'placement_phylo_placed_variants_all.tsv'),
          path(nr_matrix, stageAs: 'placement_nr_matrix.tsv'),
          path(nv_matrix, stageAs: 'placement_nv_matrix.tsv')
    val(gender)
    val(vaf_absent)
    val(vaf_present)
    val(tree_mut_pval)
    val(keep_ancestral)
    val(create_multi_tree)
    val(genotype_conv_prob)
    val(min_pval_for_true_somatic)
    val(min_variant_reads_shared)
    val(min_vaf_shared)

    output:
    tuple val(group), path("output_both_placement_*"), emit: placement_outputs

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob always resolves even when the run is skipped
    touch output_both_placement_no_results

    # ── Guard: tree file must be non-empty ────────────────────────────────────
    tree_bytes=\$(wc -c < "${tree_file}" || echo 0)
    if [ "\${tree_bytes}" -eq 0 ]; then
        echo "[both_placement] Tree file is empty (0 bytes) — skipping"
        exit 0
    fi

    NR_IN="${nr_matrix}"
    NV_IN="${nv_matrix}"
    if [ -s "${phylo_placed_variants_all}" ]; then
        tail -n +2 "${phylo_placed_variants_all}" | tr -d '\r' | awk -F'\\t' 'NF>=4 { print \$1"_"\$2"_"\$3"_"\$4 }' | sort -u > exclude_phylo_variant_ids.txt
        awk -F'\\t' 'NR==FNR{ex[\$1]=1;next} NR==1{print;next} !(\$1 in ex)' exclude_phylo_variant_ids.txt "${nr_matrix}" > nr_for_placement.tsv
        awk -F'\\t' 'NR==FNR{ex[\$1]=1;next} NR==1{print;next} !(\$1 in ex)' exclude_phylo_variant_ids.txt "${nv_matrix}" > nv_for_placement.tsv
        NR_IN="nr_for_placement.tsv"
        NV_IN="nv_for_placement.tsv"
        echo "[both_placement] Excluded \$(wc -l < exclude_phylo_variant_ids.txt) distinct variant IDs present in phylo *_placed_variants_all.tsv"
    fi

    # ── Count all variant rows in the (possibly subset) NR matrix ────────────
    n_vars=\$(tail -n +2 "\${NR_IN}" | wc -l)
    echo "[both_placement] Total variants in matrix: \${n_vars}"

    if [ "\${n_vars}" -eq 0 ]; then
        echo "[both_placement] No variants — skipping"
        exit 0
    fi

    mkdir -p output_both_placement_${group}

    cp "${tree_file}" "output_both_placement_${group}/"

    if [ -s "${phylo_placed_variants_all}" ]; then
        cp "${phylo_placed_variants_all}" "output_both_placement_${group}/"
    fi

    Rscript /usr/local/bin/rscript_minimal_variant_placement.R \\
        --donor_id                   "${group}" \\
        --input_nr                   "\${NR_IN}" \\
        --input_nv                   "\${NV_IN}" \\
        --input_tree                 "${tree_file}" \\
        --output_dir                 "output_both_placement_${group}/" \\
        --variant_type               both \\
        --tree_mut_pval              ${tree_mut_pval} \\
        --phylo_placed_variants_tsv  "${phylo_placed_variants_all}"
    """
}

// ---------------------------------------------------------------------------
// POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_{SNV,INDEL,BOTH}
// VAF + digital genotype heatmap postprocessing on placed variants.
// Input tuple: (group, placement_dirs, nr, nv, gt_files, mandatory_qc, genotype_bin).
// Emits individual res_* files (or a res_no_results_*.txt sentinel on skip).
// ---------------------------------------------------------------------------
process POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_SNV {
    tag "${group}"

    input:
    tuple val(group), path(snv_placement_dirs), path(nr_matrix), path(nv_matrix), path(gt_files), path(mandatory_qc, stageAs: 'mandatory_qc_input.tsv'), path(genotype_bin, stageAs: 'genotype_bin_input.tsv')

    output:
    tuple val(group), path("res_*_SNV.*"), emit: postprocess_snv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob (res_*_SNV.*) always resolves — even when this
    # variant type is skipped below (keeps the downstream GATHER join alive).
    # Removed on the success path before the real result files are written.
    touch res_no_results_SNV.txt

    for f in *_genotype.tsv; do
        sample="\${f%_genotype.tsv}"
        tail -n +2 "\${f}" | awk -v s="\${sample}" 'BEGIN{OFS="\\t"} {print s, \$1, \$2}'
    done > df_all_gt.tsv

    awk 'NR>1 {gsub(/"/, "", \$1); print \$1}' ${nr_matrix} > ids.txt
    awk -v OFS="\\t" -v FS="\\t" \\
        'NR==FNR {a[\$0]; next} \$2 in a {print}' \\
        ids.txt df_all_gt.tsv > df_all_gt_chosen.tsv
    echo "GT rows (filtered to NR matrix variants): \$(wc -l < df_all_gt_chosen.tsv)"

    echo "=== [snv] Postprocessing ==="
    real_dir=\$(ls -d output_snv_placement_*/ 2>/dev/null | head -1 || true)
    if [ -z "\${real_dir}" ]; then
        echo "[snv] No placement directory found — skipping"
        exit 0
    fi

    file_scheme="\${real_dir}${group}_snv_assigned_to_branches_selectedscheme.txt"
    if [ ! -s "\${file_scheme}" ]; then
        echo "[snv] Missing or empty \${file_scheme} — skipping heatmap"
        exit 0
    fi

    file_tree=\$(ls "\${real_dir}"/*_tree_with_branch_length_selectedscheme.tree 2>/dev/null | head -1 || true)
    if [ -z "\${file_tree}" ]; then
        file_tree=\$(ls "\${real_dir}"/*_tree_with_branch_length_allunfilteredvariants.tree 2>/dev/null | head -1 || true)
    fi
    if [ -z "\${file_tree}" ]; then
        echo "[snv] Missing tree Newick — skipping"
        exit 0
    fi

    n_placed=\$(awk -F'\\t' 'NR>1 { gsub(/\\r/, "", \$1); if (\$1 == "phylogeny_filtered_variant_placement") n++ } END { print n+0 }' "\${file_scheme}")
    echo "[snv] Rows with provenance phylogeny_filtered_variant_placement in \${file_scheme}: \${n_placed}"
    if [ "\${n_placed}" -eq 0 ]; then
        echo "[snv] No phylogeny_filtered_variant_placement rows — skipping heatmap"
        exit 0
    fi

    echo "[snv] Tree file for rscript_postprocess_vaf_tree_egg.R: \${file_tree}"

    Rscript /usr/local/bin/rscript_postprocess_vaf_tree_egg.R \\
        "\${file_scheme}" \\
        "${nv_matrix}" \\
        "${nr_matrix}" \\
        "\${file_tree}" \\
        df_all_gt_chosen.tsv \\
        /dev/null \\
        ${params.postprocess_heatmap_y_axis_text_size} \\
        ${params.plot_tips} \\
        "${mandatory_qc}" \\
        "${genotype_bin}"

    # Success: drop the skip sentinel and emit the result files individually.
    rm -f res_no_results_SNV.txt
    mv res_composition.pdf         res_composition_vaf_SNV.pdf
    mv res_composition_digital.pdf res_composition_digital_SNV.pdf
    mv res_figures.RDS             res_figures_SNV.RDS

    echo "[snv] Done"
    """
}

process POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_INDEL {
    tag "${group}"

    input:
    tuple val(group), path(indel_placement_dirs), path(nr_matrix), path(nv_matrix), path(gt_files), path(mandatory_qc, stageAs: 'mandatory_qc_input.tsv'), path(genotype_bin, stageAs: 'genotype_bin_input.tsv')

    output:
    tuple val(group), path("res_*_Indel.*"), emit: postprocess_indel

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob (res_*_Indel.*) always resolves — even when this
    # variant type is skipped below (keeps the downstream GATHER join alive).
    touch res_no_results_Indel.txt

    for f in *_genotype.tsv; do
        sample="\${f%_genotype.tsv}"
        tail -n +2 "\${f}" | awk -v s="\${sample}" 'BEGIN{OFS="\\t"} {print s, \$1, \$2}'
    done > df_all_gt.tsv

    awk 'NR>1 {gsub(/"/, "", \$1); print \$1}' ${nr_matrix} > ids.txt
    awk -v OFS="\\t" -v FS="\\t" \\
        'NR==FNR {a[\$0]; next} \$2 in a {print}' \\
        ids.txt df_all_gt.tsv > df_all_gt_chosen.tsv
    echo "GT rows (filtered to NR matrix variants): \$(wc -l < df_all_gt_chosen.tsv)"

    echo "=== [indel] Postprocessing ==="
    real_dir=\$(ls -d output_indel_placement_*/ 2>/dev/null | head -1 || true)
    if [ -z "\${real_dir}" ]; then
        echo "[indel] No placement directory found — skipping"
        exit 0
    fi

    file_scheme="\${real_dir}${group}_indel_assigned_to_branches_selectedscheme.txt"
    if [ ! -s "\${file_scheme}" ]; then
        echo "[indel] Missing or empty \${file_scheme} — skipping heatmap"
        exit 0
    fi

    file_tree=\$(ls "\${real_dir}"/*_tree_with_branch_length_selectedscheme.tree 2>/dev/null | head -1 || true)
    if [ -z "\${file_tree}" ]; then
        file_tree=\$(ls "\${real_dir}"/*_tree_with_branch_length_allunfilteredvariants.tree 2>/dev/null | head -1 || true)
    fi
    if [ -z "\${file_tree}" ]; then
        echo "[indel] Missing tree Newick — skipping"
        exit 0
    fi

    n_placed=\$(awk -F'\\t' 'NR>1 { gsub(/\\r/, "", \$1); if (\$1 == "phylogeny_filtered_variant_placement") n++ } END { print n+0 }' "\${file_scheme}")
    echo "[indel] Rows with provenance phylogeny_filtered_variant_placement in \${file_scheme}: \${n_placed}"
    if [ "\${n_placed}" -eq 0 ]; then
        echo "[indel] No phylogeny_filtered_variant_placement rows — skipping heatmap"
        exit 0
    fi

    echo "[indel] Tree file for rscript_postprocess_vaf_tree_egg.R: \${file_tree}"

    Rscript /usr/local/bin/rscript_postprocess_vaf_tree_egg.R \\
        "\${file_scheme}" \\
        "${nv_matrix}" \\
        "${nr_matrix}" \\
        "\${file_tree}" \\
        df_all_gt_chosen.tsv \\
        /dev/null \\
        ${params.postprocess_heatmap_y_axis_text_size} \\
        ${params.plot_tips} \\
        "${mandatory_qc}" \\
        "${genotype_bin}"

    rm -f res_no_results_Indel.txt
    mv res_composition.pdf         res_composition_vaf_Indel.pdf
    mv res_composition_digital.pdf res_composition_digital_Indel.pdf
    mv res_figures.RDS             res_figures_Indel.RDS

    echo "[indel] Done"
    """
}

process POSTPROCESS_SEQUOIA_DRAWVAFHEAT_TREE_BOTH {
    tag "${group}"

    input:
    tuple val(group), path(both_placement_dirs), path(nr_matrix), path(nv_matrix), path(gt_files), path(mandatory_qc, stageAs: 'mandatory_qc_input.tsv'), path(genotype_bin, stageAs: 'genotype_bin_input.tsv')

    output:
    tuple val(group), path("res_*_BOTH.*"), emit: postprocess_both

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Sentinel so the output glob (res_*_BOTH.*) always resolves — even when this
    # variant type is skipped below (keeps the downstream GATHER join alive).
    touch res_no_results_BOTH.txt

    for f in *_genotype.tsv; do
        sample="\${f%_genotype.tsv}"
        tail -n +2 "\${f}" | awk -v s="\${sample}" 'BEGIN{OFS="\\t"} {print s, \$1, \$2}'
    done > df_all_gt.tsv

    awk 'NR>1 {gsub(/"/, "", \$1); print \$1}' ${nr_matrix} > ids.txt
    awk -v OFS="\\t" -v FS="\\t" \\
        'NR==FNR {a[\$0]; next} \$2 in a {print}' \\
        ids.txt df_all_gt.tsv > df_all_gt_chosen.tsv
    echo "GT rows (filtered to NR matrix variants): \$(wc -l < df_all_gt_chosen.tsv)"

    echo "=== [both] Postprocessing ==="
    real_dir=\$(ls -d output_both_placement_*/ 2>/dev/null | head -1 || true)
    if [ -z "\${real_dir}" ]; then
        echo "[both] No placement directory found — skipping"
        exit 0
    fi

    file_scheme="\${real_dir}${group}_both_assigned_to_branches_selectedscheme.txt"
    if [ ! -s "\${file_scheme}" ]; then
        echo "[both] Missing or empty \${file_scheme} — skipping heatmap"
        exit 0
    fi

    file_tree=\$(ls "\${real_dir}"/*_tree_with_branch_length_selectedscheme.tree 2>/dev/null | head -1 || true)
    if [ -z "\${file_tree}" ]; then
        file_tree=\$(ls "\${real_dir}"/*_tree_with_branch_length_allunfilteredvariants.tree 2>/dev/null | head -1 || true)
    fi
    if [ -z "\${file_tree}" ]; then
        echo "[both] Missing tree Newick — skipping"
        exit 0
    fi

    n_placed=\$(awk -F'\\t' 'NR>1 { gsub(/\\r/, "", \$1); if (\$1 == "phylogeny_filtered_variant_placement") n++ } END { print n+0 }' "\${file_scheme}")
    echo "[both] Rows with provenance phylogeny_filtered_variant_placement in \${file_scheme}: \${n_placed}"
    if [ "\${n_placed}" -eq 0 ]; then
        echo "[both] No phylogeny_filtered_variant_placement rows — skipping heatmap"
        exit 0
    fi

    echo "[both] Tree file for rscript_postprocess_vaf_tree_egg.R: \${file_tree}"

    Rscript /usr/local/bin/rscript_postprocess_vaf_tree_egg.R \\
        "\${file_scheme}" \\
        "${nv_matrix}" \\
        "${nr_matrix}" \\
        "\${file_tree}" \\
        df_all_gt_chosen.tsv \\
        /dev/null \\
        ${params.postprocess_heatmap_y_axis_text_size} \\
        ${params.plot_tips} \\
        "${mandatory_qc}" \\
        "${genotype_bin}"

    rm -f res_no_results_BOTH.txt
    mv res_composition.pdf         res_composition_vaf_BOTH.pdf
    mv res_composition_digital.pdf res_composition_digital_BOTH.pdf
    mv res_figures.RDS             res_figures_BOTH.RDS

    echo "[both] Done"
    """
}

// ---------------------------------------------------------------------------
// GATHER_LINEAGE_POSTPROCESS_ARTIFACTS
// Collect per-group postprocess figures + the placed-variant tables/trees they
// were built from into one lineage_artifacts_${group}/ folder. The phylogeny
// output dirs are also staged so their tree PDFs (tree_with_branch_length*.pdf)
// are collected and rendered to PNG here — giving the bundle (and MultiQC) tree
// figures even when the placement/heatmap path is skipped.
// Input tuple: (group, postprocess_files, placement_dirs, phylogeny_dirs).
// ---------------------------------------------------------------------------
process GATHER_LINEAGE_POSTPROCESS_ARTIFACTS {
    tag "${group}"

    input:
    tuple val(group), path(postprocess_files), path(placement_dirs), path(phylogeny_dirs)

    output:
    tuple val(group), path("lineage_artifacts_${group}"), emit: bundle

    script:
    """
    set -euo pipefail
    dest="lineage_artifacts_${group}"
    mkdir -p "\$dest/figures_pdf" "\$dest/figures_png" "\$dest/placed_variants" "\$dest/trees"

    # 1) PDFs created by the postprocess. Prune \$dest so we never re-collect our own output.
    #    -L: inputs are staged as SYMLINKS; follow them into dirs / resolve file links.
    find -L . -path "./\$dest" -prune -o -name '*.pdf' -print 2>/dev/null | while read -r f; do
        cp -f "\$f" "\$dest/figures_pdf/" 2>/dev/null || true
    done

    # 2) PNGs: copy any that already exist, then render one PNG per gathered PDF.
    find -L . -path "./\$dest" -prune -o -name '*.png' -print 2>/dev/null | while read -r f; do
        cp -f "\$f" "\$dest/figures_png/" 2>/dev/null || true
    done
    for pdf in "\$dest"/figures_pdf/*.pdf; do
        [ -e "\$pdf" ] || continue
        base=\$(basename "\$pdf" .pdf)
        pdftoppm -png -r 150 "\$pdf" "\$dest/figures_png/\${base}" >/dev/null 2>&1 || true
    done

    # 3) Placed-variant tables the heatmaps were built from.
    find -L . -path "./\$dest" -prune -o -name '*_assigned_to_branches_selectedscheme.txt' -print 2>/dev/null | while read -r f; do
        cp -f "\$f" "\$dest/placed_variants/" 2>/dev/null || true
    done

    # 4) Trees actually used for the tree panel, one per placement / phylogeny dir.
    for pdir in output_*_placement_* output_*_unfiltered; do
        [ -d "\$pdir" ] || continue
        t=\$(ls "\$pdir"/*_tree_with_branch_length_selectedscheme.tree 2>/dev/null | head -1 || true)
        if [ -z "\$t" ]; then
            t=\$(ls "\$pdir"/*_tree_with_branch_length_allunfilteredvariants.tree 2>/dev/null | head -1 || true)
        fi
        if [ -z "\$t" ]; then
            t=\$(ls "\$pdir"/*_tree_with_branch_length.tree 2>/dev/null | head -1 || true)
        fi
        [ -n "\$t" ] && cp -f "\$t" "\$dest/trees/" 2>/dev/null || true
    done

    n_pdf=\$(ls -1 "\$dest/figures_pdf" 2>/dev/null | wc -l)
    n_png=\$(ls -1 "\$dest/figures_png" 2>/dev/null | wc -l)
    n_var=\$(ls -1 "\$dest/placed_variants" 2>/dev/null | wc -l)
    n_tree=\$(ls -1 "\$dest/trees" 2>/dev/null | wc -l)
    echo "[gather] ${group}: \${n_pdf} pdf, \${n_png} png, \${n_var} placed-variant tables, \${n_tree} trees"
    find "\$dest" -type f | sort
    if [ "\$((n_pdf + n_png + n_var + n_tree))" -eq 0 ]; then
        echo "[gather] WARN: nothing collected for ${group}" >&2
    fi
    """
}

// ============================================================================
// MUTSIGNATURES BRANCH
// ============================================================================
// Ported (basej style) from $REF:modules/bioskryb/<dir>/main.nf. All publishDir,
// container, val(publish_dir)/val(enable_publish) inputs and params.timestamp
// usage have been stripped. Containers/resources live in nextflow.config.
// ============================================================================

// ---------------------------------------------------------------------------
// EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE_FROM_BINARY_MATRIX
// Generate per-sample mutsig somatic variant tables from a binary prevalence
// matrix. Row keys must follow {chrom}_{pos}_{ref}_{alt} (pos = rightmost numeric
// token). Only rows with value == 1 are written for each sample. Output schema
// matches the VCF-mode exporter's schema and feeds SigProfilerMatrixGenerator.
// ---------------------------------------------------------------------------
process EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE_FROM_BINARY_MATRIX {
    tag "${cohort_id}"

    input:
    path(binary_matrix)
    val(cohort_id)
    val(genome_build)

    output:
    path("*_mutsig_somatic_variant_table.tsv"), emit: per_sample_tsvs

    script:
    """
    set -euo pipefail
    python3 << 'PY'
import os

HEADER = 'Project\\tSample\\tID\\tGenome\\tmut_type\\tchrom\\tpos_start\\tpos_end\\tref\\talt\\tType'
PROJECT  = '${cohort_id}'
GENOME   = '${genome_build}'

def parse_variant_id(vid):
    # Format: {chrom}_{pos}_{ref}_{alt}
    # pos is the rightmost numeric token; everything before it is chrom.
    parts = vid.split('_')
    # Walk from right: last=alt, second-to-last=ref, third-to-last=pos (numeric)
    alt = parts[-1]
    ref = parts[-2]
    pos = parts[-3]
    chrom = '_'.join(parts[:-3])
    return chrom, pos, ref, alt

def mut_type_and_end(ref, alt, pos):
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP', pos
    elif len(ref) == 2 and len(alt) == 2:
        return 'DBS', pos
    else:
        return 'INDEL', str(int(pos) + len(ref) - 1)

with open('${binary_matrix}') as fh:
    lines = [l.rstrip('\\n') for l in fh if l.strip()]

header_parts = lines[0].split('\\t')
samples = header_parts[1:]  # column names

# Open one output file per sample
out_handles = {}
for s in samples:
    fname = f'{s}_mutsig_somatic_variant_table.tsv'
    fh = open(fname, 'w')
    fh.write(HEADER + '\\n')
    out_handles[s] = fh

counts = {s: 0 for s in samples}

for line in lines[1:]:
    parts = line.split('\\t')
    vid  = parts[0]
    bits = parts[1:]
    chrom, pos, ref, alt = parse_variant_id(vid)
    mtype, pos_end = mut_type_and_end(ref, alt, pos)
    row = f'{PROJECT}\\t{{sample}}\\t.\\t{GENOME}\\t{mtype}\\t{chrom}\\t{pos}\\t{pos_end}\\t{ref}\\t{alt}\\tSOMATIC'
    for i, s in enumerate(samples):
        if i < len(bits) and bits[i] == '1':
            out_handles[s].write(row.format(sample=s) + '\\n')
            counts[s] += 1

for fh in out_handles.values():
    fh.close()

total = sum(counts.values())
nonempty = sum(1 for v in counts.values() if v > 0)
print(
    f'[export_mutsig_somatic_variant_table_from_binary_matrix] '
    f'{len(samples)} samples, {nonempty} with >=1 variant, {total} total variant-sample rows written',
    flush=True
)
PY
    """
}

// ---------------------------------------------------------------------------
// SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES
// Collect per-sample variant TSVs into mat_files/ (.txt copies; SigProfiler rejects
// .tsv) and run SigProfilerMatrixGeneratorFunc (SBS96, DBS78, ID83).
// ---------------------------------------------------------------------------
process SIGPROFILER_MATRIX_GENERATOR_FROM_MUTSIG_SOMATIC_VARIANT_TABLES {
    tag "${cohort_id}"

    input:
    path tsv_inputs
    val(cohort_id)
    val(sig_genome_build)
    val(sig_exome)

    output:
    path('matrix_generator_out'), emit: matrix_output

    script:
    def exomePy = sig_exome.toString() == 'true' ? 'True' : 'False'
    """
    set -euo pipefail
    python3 <<'PY'
import glob
import os
import shutil

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

project = "${cohort_id}"
reference_genome = "${sig_genome_build}"
exome = ${exomePy}
tsv_glob = "*_mutsig_somatic_variant_table.tsv"
out_dir_name = "matrix_generator_out"

os.makedirs("mat_files", exist_ok=True)
paths = sorted(glob.glob(tsv_glob))
if not paths:
    raise SystemExit("No TSV files matched: " + tsv_glob)
for p in paths:
    base = os.path.basename(p)
    if base.endswith(".tsv"):
        dst = os.path.join("mat_files", base[:-4] + ".txt")
    else:
        dst = os.path.join("mat_files", base + ".txt")
    shutil.copy2(p, dst)

mat_dir = os.path.abspath("mat_files") + "/"
out_dir = os.path.abspath(out_dir_name) + "/"
matGen.SigProfilerMatrixGeneratorFunc(
    project,
    reference_genome,
    mat_dir,
    exome=exome,
    bed_file=None,
    chrom_based=False,
    plot=False,
    tsb_stat=False,
    seqInfo=False,
    output_directory=out_dir,
)
PY
    """
}

// ---------------------------------------------------------------------------
// SIGPROFILER_ASSIGNMENT_FROM_MATRIX
// Cohort-level COSMIC signature assignment from a SigProfilerMatrixGenerator matrix.
// Supports SBS96 / DBS78 / ID83 via the mutation_type input ("SBS"/"DBS"/"ID").
// ---------------------------------------------------------------------------
process SIGPROFILER_ASSIGNMENT_FROM_MATRIX {
    tag "${cohort_id}"

    input:
    path(matrix_generator_out)
    val(mutation_type)
    val(cohort_id)
    val(genome_build)
    val(context_type)
    val(exome)
    val(export_probabilities)
    val(export_probabilities_per_mutation)

    output:
    path("${cohort_id}_${mutation_type}_sig_assignment"), emit: assignment_cohort_dir

    script:
    def exomePy = exome.toString() == 'true' ? 'True' : 'False'
    def exportProbPy = export_probabilities.toString() == 'true' ? 'True' : 'False'
    def exportPerMutPy = export_probabilities_per_mutation.toString() == 'true' ? 'True' : 'False'
    // Must match MutationMatrixGenerator output suffix (see SigProfilerMatrixGeneratorFunc.py).
    def matrixScopeSuffix = exome.toString() == 'true' ? '.exome' : '.all'
    """
    set -euo pipefail
    python3 <<'PY'
import os
from SigProfilerAssignment import Analyzer as A

cohort = "${cohort_id}"
mutation_type = "${mutation_type}"
scope_suffix = "${matrixScopeSuffix}"

# Determine subdirectory, matrix filename base, and context_type based on mutation_type
if mutation_type == "SBS":
    subdir = "SBS"
    ctx = "${context_type}"
    base = cohort + ".SBS" + ctx + scope_suffix
    fit_context_type = "${context_type}"
    collapse_to_SBS96 = True
elif mutation_type == "DBS":
    subdir = "DBS"
    base = cohort + ".DBS78" + scope_suffix
    fit_context_type = "DBS78"
    collapse_to_SBS96 = False
elif mutation_type == "ID":
    subdir = "ID"
    base = cohort + ".ID83" + scope_suffix
    fit_context_type = "ID83"
    collapse_to_SBS96 = False
else:
    raise SystemExit("Unknown mutation_type: " + mutation_type)

samples = os.path.abspath(os.path.join("matrix_generator_out", subdir, base))
out = os.path.abspath(cohort + "_" + mutation_type + "_sig_assignment")
if not os.path.isfile(samples):
    raise SystemExit("Missing " + mutation_type + " matrix: " + samples)
os.makedirs(out, exist_ok=True)
A.cosmic_fit(
    samples=samples,
    output=out,
    input_type="matrix",
    context_type=fit_context_type,
    collapse_to_SBS96=collapse_to_SBS96,
    genome_build="${genome_build}",
    exome=${exomePy},
    make_plots=False,
    export_probabilities=${exportProbPy},
    export_probabilities_per_mutation=${exportPerMutPy},
    verbose=False,
)
PY
    """
}

// ---------------------------------------------------------------------------
// MERGE_SIGNATURE_ACTIVITIES
// Merge per-sample Assignment_Solution_Activities.txt into a cohort-level
// SamplesXSignatures matrix. Parameterized by mutation_type (SBS/DBS/ID).
// ---------------------------------------------------------------------------
process MERGE_SIGNATURE_ACTIVITIES {

    input:
    path(assignment_dirs)   // collected list of per-sample SigProfilerAssignment output directories
    val(mutation_type)      // one of "SBS", "DBS", "ID" — used to differentiate output filenames

    output:
    path("merged_signature_activities_${mutation_type}.txt"), emit: merged_activities

    script:
    """
    python /scripts/merge_signature_activities.py \\
        --input_dirs ${assignment_dirs} \\
        --output     merged_signature_activities_${mutation_type}.txt
    """
}

// ---------------------------------------------------------------------------
// COMPUTE_MUTSIG_COVERAGE_FROM_VARIANT_TABLE
// Per-sample mutational signature coverage from a variant table TSV.
//   InputVariants                    = data rows in the per-sample variant table
//   NumberVariantsAssignedMutSig     = sum of activity counts (SBS + DBS + ID dirs)
//   ProportionVariantsAssignedMutSig = assigned / input
// ---------------------------------------------------------------------------
process COMPUTE_MUTSIG_COVERAGE_FROM_VARIANT_TABLE {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(variant_table), path(assignment_dirs)

    output:
    tuple val(sample_name), path("mutsig_coverage_${sample_name}.tsv"), emit: coverage_tsv

    script:
    """
    set -euo pipefail

    # Total variants for this sample (header row excluded)
    input_vars=\$(tail -n +2 "${variant_table}" | wc -l)

    # Sum activity counts across all SigProfilerAssignment directories (SBS + DBS + ID).
    total_assigned=0
    for dir in ${assignment_dirs}; do
        acts_file="\${dir}/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
        if [ -f "\$acts_file" ]; then
            sample_sum=\$(awk -v sid='${sample_name}' 'BEGIN{FS=OFS="\\t"} NR==1{next} \$1==sid {s=0; for(i=2;i<=NF;i++) s+=\$i; print s; exit}' "\$acts_file")
            if [ -n "\$sample_sum" ]; then
                total_assigned=\$((total_assigned + \${sample_sum%.*}))
            fi
        fi
    done

    prop=\$(awk -v inp="\${input_vars}" -v asgn="\${total_assigned}" \\
        'BEGIN { print (inp > 0) ? asgn/inp : 0 }')

    printf "SampleID\\tInputVariants\\tNumberVariantsAssignedMutSig\\tProportionVariantsAssignedMutSig\\n" \\
        > mutsig_coverage_${sample_name}.tsv
    printf "%s\\t%s\\t%s\\t%s\\n" "${sample_name}" "\${input_vars}" "\${total_assigned}" "\${prop}" \\
        >> mutsig_coverage_${sample_name}.tsv
    """
}

// ---------------------------------------------------------------------------
// PLOT_MATRIX_SIGNATURE_BARGRAPHS
// Signature activity bar graphs for matrix-based SigProfilerAssignment
// (merged_signature_activities_{mutation_type}.txt). Supports SBS/DBS/ID.
// Emits combined PNG + PDF (title includes cohort % assigned) and interactive HTML.
// ---------------------------------------------------------------------------
process PLOT_MATRIX_SIGNATURE_BARGRAPHS {
    tag 'matrix_sig_bars'

    input:
    tuple path(merged_signature_activities), path(mutsig_coverage_tsvs)
    val(mutation_type)

    output:
    path("signature_bargraphs_${mutation_type}_combined.png"),     emit: bargraph_png
    path("signature_bargraphs_${mutation_type}_combined.pdf"),     emit: bargraph_pdf
    path("signature_bargraphs_${mutation_type}_interactive.html"), emit: bargraph_html

    script:
    """
    set -euo pipefail
    shopt -s nullglob
    cov_files=(mutsig_coverage_*.tsv)
    if [ \${#cov_files[@]} -eq 0 ]; then
      echo "PLOT_MATRIX_SIGNATURE_BARGRAPHS: no mutsig_coverage_*.tsv in task directory" >&2
      exit 1
    fi
    head -1 "\${cov_files[0]}" > mutsig_coverage_merged.tsv
    for f in "\${cov_files[@]}"; do tail -n +2 "\$f" >> mutsig_coverage_merged.tsv; done

    Rscript /usr/local/bin/plot_matrix_signature_bargraphs.R \\
      ${merged_signature_activities} \\
      mutsig_coverage_merged.tsv

    mv signature_bargraphs_combined.png signature_bargraphs_${mutation_type}_combined.png
    mv signature_bargraphs_combined.pdf signature_bargraphs_${mutation_type}_combined.pdf

    Rscript /usr/local/bin/generate_interactive_matrix_bargraph.R \\
      ${merged_signature_activities}

    mv signature_bargraphs_interactive.html signature_bargraphs_${mutation_type}_interactive.html
    """
}

// ---------------------------------------------------------------------------
// MUSICAL_SIGNATURE_ANALYSIS
// MuSiCaL (Park lab) COSMIC v3.2 SBS refitting on the cohort SBS96 matrix
// (likelihood-bidirectional / sparse NNLS). Supports WES / WGS via exome flag.
// ---------------------------------------------------------------------------
process MUSICAL_SIGNATURE_ANALYSIS {
    tag "${cohort_id}"

    input:
    path(matrix_generator_out)
    val(cohort_id)
    val(genome_build)
    val(exome)

    output:
    path("${cohort_id}_musical_activities.tsv"),            emit: musical_activities
    path("${cohort_id}_musical_cosine_similarities.tsv"),   emit: musical_cosine_sims

    script:
    def exomeStr         = exome.toString() == 'true' ? 'true' : 'false'
    def matrixScopeSuffix = exome.toString() == 'true' ? 'exome' : 'all'
    """
    set -euo pipefail
    matrix_file="${matrix_generator_out}/SBS/${cohort_id}.SBS96.${matrixScopeSuffix}"
    if [ ! -f "\${matrix_file}" ]; then
        echo "MUSICAL_SIGNATURE_ANALYSIS: SBS96 matrix not found: \${matrix_file}" >&2
        exit 1
    fi
    python3 /usr/local/bin/run_musical_analysis.py \\
        "\${matrix_file}" \\
        "${cohort_id}" \\
        "${genome_build}" \\
        "${exomeStr}" \\
        "${cohort_id}"
    """
}

// ---------------------------------------------------------------------------
// SIGDYN_SIGNATURE_ANALYSIS
// SigDyn (Goncalves lab) dynamics analysis via R/MutationalPatterns strict
// best-subset NNLS COSMIC v3 SBS fitting + per-sample cosine similarity +
// signature dynamics statistics. Takes the matrix-generator dir.
// ---------------------------------------------------------------------------
process SIGDYN_SIGNATURE_ANALYSIS {
    tag "${cohort_id}"

    input:
    path(matrix_generator_out)
    val(cohort_id)
    val(genome_build)
    val(exome)

    output:
    path("${cohort_id}_sigdyn_activities.tsv"),             emit: sigdyn_activities
    path("${cohort_id}_sigdyn_cosine_similarities.tsv"),    emit: sigdyn_cosine_sims
    path("${cohort_id}_sigdyn_signature_stats.tsv"),        emit: sigdyn_signature_stats
    path("${cohort_id}_sigdyn_activity_heatmap.png"),       emit: sigdyn_heatmap,  optional: true
    path("${cohort_id}_sigdyn_signature_dynamics_cv.png"),  emit: sigdyn_cv_plot,  optional: true

    script:
    def exomeStr          = exome.toString() == 'true' ? 'true' : 'false'
    def matrixScopeSuffix = exome.toString() == 'true' ? 'exome' : 'all'
    """
    set -euo pipefail
    matrix_file="${matrix_generator_out}/SBS/${cohort_id}.SBS96.${matrixScopeSuffix}"
    if [ ! -f "\${matrix_file}" ]; then
        echo "SIGDYN_SIGNATURE_ANALYSIS: SBS96 matrix not found: \${matrix_file}" >&2
        exit 1
    fi
    Rscript /usr/local/bin/run_sigdyn_analysis.R \\
        "\${matrix_file}" \\
        "${cohort_id}" \\
        "${genome_build}" \\
        "${exomeStr}" \\
        "${cohort_id}"
    """
}

// ---------------------------------------------------------------------------
// MUSICATK_SIGNATURE_ANALYSIS
// musicatk (Campbell lab) de novo NMF/LDA signature discovery from per-sample
// mutsig variant TSVs. Requires R >= 4.4.0 (disabled by default).
//   musicatk_k_denovo: number of de novo signatures to extract.
// ---------------------------------------------------------------------------
process MUSICATK_SIGNATURE_ANALYSIS {
    tag "${cohort_id}"

    input:
    path(per_sample_tsvs)    // collected list of *_mutsig_somatic_variant_table.tsv
    val(cohort_id)
    val(genome_build)
    val(exome)
    val(k_denovo)

    output:
    path("${cohort_id}_musicatk_denovo_activities.tsv"),    emit: musicatk_activities
    path("${cohort_id}_musicatk_denovo_signatures.tsv"),    emit: musicatk_signatures
    path("${cohort_id}_musicatk_denovo_exposures.png"),     emit: musicatk_exposures_png,   optional: true
    path("${cohort_id}_musicatk_denovo_signatures.png"),    emit: musicatk_signatures_png,  optional: true

    script:
    def exomeStr = exome.toString() == 'true' ? 'true' : 'false'
    """
    set -euo pipefail
    Rscript /usr/local/bin/run_musicatk_analysis.R \\
        "${per_sample_tsvs}" \\
        "${cohort_id}" \\
        "${genome_build}" \\
        "${exomeStr}" \\
        "${k_denovo}" \\
        "${cohort_id}"
    """
}

// ---------------------------------------------------------------------------
// PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES
// Comparative faceted bar chart (SigProfilerAssignment, MuSiCaL, SigDyn) sharing
// x-axis sample order + signature palette; plus concordance figures. Requires
// activity matrices from all three methods.
// ---------------------------------------------------------------------------
process PLOT_COMPARATIVE_SIGNATURE_ACTIVITIES {
    tag "${cohort_id}"

    input:
    path(spa_activities)      // merged_signature_activities_SBS.txt
    path(musical_activities)  // {cohort_id}_musical_activities.tsv
    path(sigdyn_activities)   // {cohort_id}_sigdyn_activities.tsv
    val(cohort_id)

    output:
    path("${cohort_id}_comparative_signature_activities.png"),  emit: comparative_png
    path("${cohort_id}_comparative_signature_activities.pdf"),  emit: comparative_pdf
    path("${cohort_id}_concordance_detection.png"),             emit: concordance_detection
    path("${cohort_id}_concordance_sample.png"),                emit: concordance_sample
    path("${cohort_id}_concordance_signature_corr.png"),        emit: concordance_sig_corr

    script:
    """
    set -euo pipefail
    Rscript /usr/local/bin/plot_comparative_signature_activities.R \\
        "${spa_activities}" \\
        "${musical_activities}" \\
        "${sigdyn_activities}" \\
        "${cohort_id}" \\
        "${cohort_id}"
    """
}

// ---------------------------------------------------------------------------
// GATHER_MUTSIG_ARTIFACTS
// Collect, per cohort, the key mutational-signature deliverables into one
// mutsig_artifacts_${cohort_id}/ folder (activity_matrices/, bargraph_figures/,
// comparative_figures/). Any category may be empty; the script tolerates it.
// ---------------------------------------------------------------------------
process GATHER_MUTSIG_ARTIFACTS {
    tag "${cohort_id}"

    input:
    val(cohort_id)
    path(activity_matrices,   stageAs: "in_activity/*")
    path(bargraph_figures,    stageAs: "in_bargraphs/*")
    path(comparative_figures, stageAs: "in_comparative/*")

    output:
    tuple val(cohort_id), path("mutsig_artifacts_${cohort_id}"), emit: bundle

    script:
    """
    set -euo pipefail
    dest="mutsig_artifacts_${cohort_id}"
    mkdir -p "\$dest/activity_matrices" "\$dest/bargraph_figures" "\$dest/comparative_figures"

    # Merged activity matrices for all methods (SPA SBS/DBS/ID, MuSiCaL, SigDyn).
    [ -d in_activity ]    && cp -f in_activity/*    "\$dest/activity_matrices/"   2>/dev/null || true
    # Signature bargraph figures (PNG + PDF) for every mutation type that ran (SBS/DBS/ID).
    [ -d in_bargraphs ]   && cp -f in_bargraphs/*   "\$dest/bargraph_figures/"    2>/dev/null || true
    # Comparative + concordance figures.
    [ -d in_comparative ] && cp -f in_comparative/* "\$dest/comparative_figures/" 2>/dev/null || true

    n_act=\$(ls -1 "\$dest/activity_matrices"   2>/dev/null | wc -l)
    n_bar=\$(ls -1 "\$dest/bargraph_figures"    2>/dev/null | wc -l)
    n_cmp=\$(ls -1 "\$dest/comparative_figures" 2>/dev/null | wc -l)
    echo "[gather] ${cohort_id}: \${n_act} activity matrices, \${n_bar} bargraph figures, \${n_cmp} comparative figures"
    find "\$dest" -type f | sort
    if [ "\$((n_act + n_bar + n_cmp))" -eq 0 ]; then
        echo "[gather] WARN: nothing collected for ${cohort_id}" >&2
    fi
    """
}
