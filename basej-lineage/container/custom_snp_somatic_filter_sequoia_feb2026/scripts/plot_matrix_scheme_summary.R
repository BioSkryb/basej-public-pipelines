#!/usr/bin/env Rscript
## Render a PDF with:
##   Page 1 — table of SNV / INDEL / Total counts per scheme (cohort level)
##   Page 2 — boxplots of per-sample SNV and INDEL recovery per scheme (log10 scale)
##   Page 3 — table of per-scheme median ± IQR per sample
## Called from PLOT_MATRIX_SCHEME_SUMMARY.
## Usage: Rscript plot_matrix_scheme_summary.R <group_id> <scheme_tsv> <per_sample_tsv> [upstream_per_sample_tsv]

library(gridExtra)
library(grid)
library(ggplot2)
library(dplyr)

args     <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: plot_matrix_scheme_summary.R <group_id> <scheme_tsv> [per_sample_tsv] [upstream_per_sample_tsv]")
group_id    <- args[1]
tsv_file    <- args[2]
ps_file     <- if (length(args) >= 3) args[3] else NULL
us_file     <- if (length(args) >= 4) args[4] else NULL
pdf_file    <- paste0("matrix_scheme_summary_", group_id, ".pdf")

# ── Scheme display names (page 1 table + pages 2–3 boxplots / median table) ───
# Order on figures follows scheme_order below (mirrors the pipeline filtering order):
#   Input → VEP germline → Bulk germline → 1st-pass statistical →
#   HQ statistical → +Pileup QC → +Depth → +Phylogeny coverage
# The last four are the VAF-split cascade schemes emitted by
# CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF (build_vaf_split_variants_hexbin.R).
# unfiltered is omitted from page 1.
scheme_label_map <- c(
    pre_bulk   = "Input variants",
    post_vep   = "Post VEP germline filter",
    post_bulk  = "Post bulk germline filter",
    post_binom = "Post 1st-pass statistical filter",
    HQRoundStatisticalFiltered                                           = "HQ statistical filter",
    HQRoundStatisticalFilteredPlusQCFiltered                             = "HQ statistical +\npileup QC",
    HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered            = "HQ statistical +\nQC + depth",
    HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny = "HQ statistical +\nQC + depth + phylogeny coverage"
)

# ── Page 1: cohort-level summary table ───────────────────────────────────────
df <- read.delim(tsv_file, stringsAsFactors=FALSE, check.names=FALSE)
df <- df[df$scheme != "unfiltered", ]   # MatrixInput excluded from display
df$Total  <- df$NumberOfSNVs + df$NumberOfIndels
df$PctSNV <- ifelse(df$Total > 0,
                    sprintf("%.1f%%", df$NumberOfSNVs / df$Total * 100),
                    "--")
df$scheme_display <- ifelse(df$scheme %in% names(scheme_label_map),
                            scheme_label_map[df$scheme], df$scheme)
# Enforce canonical display order
df$scheme <- factor(df$scheme, levels=intersect(names(scheme_label_map), df$scheme))
df        <- df[order(df$scheme), ]

display <- data.frame(
    Scheme  = df$scheme_display,
    SNVs    = formatC(df$NumberOfSNVs,   format="d", big.mark=","),
    INDELs  = formatC(df$NumberOfIndels, format="d", big.mark=","),
    Total   = formatC(df$Total,          format="d", big.mark=","),
    `% SNV` = df$PctSNV,
    check.names      = FALSE,
    stringsAsFactors = FALSE
)

n         <- nrow(display)
row_fills <- rep(c("white", "#EEF3FB"), length.out=n)

tt <- ttheme_default(
    core    = list(fg_params=list(cex=0.85),
                   bg_params=list(fill=row_fills, col="gray85")),
    colhead = list(fg_params=list(cex=0.90, fontface="bold", col="white"),
                   bg_params=list(fill="#2166AC", col="#2166AC"))
)

pdf(pdf_file, width=11, height=5.5)
grid.newpage()

pushViewport(viewport(x=0.5, y=0.92, width=1, height=0.12, just=c("centre","centre")))
grid.text(
    paste0("NR/NV Matrix Filtering Scheme Summary - Group: ", group_id),
    gp=gpar(fontsize=14, fontface="bold", col="#2166AC")
)
popViewport()

pushViewport(viewport(x=0.5, y=0.83, width=1, height=0.08, just=c("centre","centre")))
grid.text(
    "SNV vs INDEL counts for each NR/NV matrix filtering scheme",
    gp=gpar(fontsize=9, col="gray40")
)
popViewport()

pushViewport(viewport(x=0.5, y=0.40, width=0.94, height=0.68, just=c("centre","centre")))
grid.draw(tableGrob(display, rows=NULL, theme=tt))
popViewport()

# ── Page 2: per-sample boxplots ───────────────────────────────────────────────
if (!is.null(ps_file) && file.exists(ps_file) && file.info(ps_file)$size > 0) {
    ps <- read.delim(ps_file, stringsAsFactors=FALSE, check.names=FALSE)

    # Prepend upstream filter stages if available
    if (!is.null(us_file) && file.exists(us_file) && file.info(us_file)$size > 0) {
        us <- read.delim(us_file, stringsAsFactors=FALSE, check.names=FALSE)
        ps <- rbind(us, ps)
    }

    scheme_order <- c("pre_bulk","post_vep","post_bulk","post_binom",
                      "HQRoundStatisticalFiltered",
                      "HQRoundStatisticalFilteredPlusQCFiltered",
                      "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFiltered",
                      "HQRoundStatisticalFilteredPlusQCFilteredPlusDepthFilteredForPhylogeny")

    label_map <- scheme_label_map

    # Rows with scheme not in label_map become NA in factor() → spurious "NA" on x-axis.
    # matrix_per_sample_summary always includes "unfiltered"; page 1 drops it — same here.
    ps <- ps[!is.na(ps$scheme) & nzchar(as.character(ps$scheme)) & ps$scheme %in% names(label_map), , drop=FALSE]
    if (nrow(ps) == 0) {
        warning("[plot_matrix_scheme_summary] No per-sample rows left after scheme filter; skipping pages 2–3.")
    } else {

    ps$scheme <- factor(ps$scheme, levels=intersect(scheme_order, unique(ps$scheme)))
    ps$scheme_label <- factor(
        label_map[as.character(ps$scheme)],
        levels = label_map[levels(ps$scheme)]
    )

    scheme_pal <- c("#7B3F00","#E8601C","#F6C141","#4DAF4A",  # upstream: pre_bulk, post_vep, post_bulk, post_binom
                    "#1082A2","#6A3D9A","#A0CC2C","#DD14D3")   # cascade: HQStat, +Depth, +QC, +ForPhylogeny
    fill_cols  <- scheme_pal[seq_len(nlevels(ps$scheme_label))]
    n_samples  <- length(unique(ps$sample))

    make_boxplot <- function(data, y_var, y_label, title_label) {
        ggplot(data, aes(x=scheme_label, y=.data[[y_var]], fill=scheme_label)) +
            geom_boxplot(outlier.shape=NA, alpha=0.75, width=0.55,
                         color="gray30", linewidth=0.4) +
            geom_jitter(width=0.15, size=1.6, alpha=0.65, color="#2C3E50") +
            scale_fill_manual(values=fill_cols) +
            scale_y_log10(labels=scales::label_comma()) +
            labs(title=title_label,
                 subtitle=paste0("Group: ", group_id, "  |  n=", n_samples, " samples"),
                 x=NULL, y=paste0(y_label, " (log10 scale)")) +
            theme_bw(base_size=11) +
            theme(
                legend.position    = "none",
                plot.title         = element_text(face="bold", color="#2166AC", size=13),
                plot.subtitle      = element_text(color="gray40", size=9),
                axis.text.x        = element_text(size=8, lineheight=0.85),
                axis.title.y       = element_text(size=10),
                panel.grid.major.x = element_blank(),
                panel.border       = element_rect(color="gray70")
            )
    }

    p_snv   <- make_boxplot(ps, "NumberOfSNVs",   "SNVs per sample",
                            "Per-sample SNV recovery by filtering scheme")
    p_indel <- make_boxplot(ps, "NumberOfIndels", "INDELs per sample",
                            "Per-sample INDEL recovery by filtering scheme")

    grid.newpage()
    vp_top <- viewport(x=0.5, y=0.75, width=1, height=0.5)
    vp_bot <- viewport(x=0.5, y=0.25, width=1, height=0.5)
    print(p_snv,   vp=vp_top)
    print(p_indel, vp=vp_bot)

    # ── Page 3: per-scheme median ± IQR table ─────────────────────────────────
    mean_tbl <- do.call(rbind, lapply(levels(ps$scheme), function(s) {
        sub  <- ps[ps$scheme == s, ]
        lbl  <- label_map[s]
        data.frame(
            Stage             = lbl,
            `Median SNVs`     = sprintf("%.1f", median(sub$NumberOfSNVs,                       na.rm=TRUE)),
            `IQR SNVs`        = sprintf("%.1f", IQR(sub$NumberOfSNVs,                          na.rm=TRUE)),
            `Median INDELs`   = sprintf("%.1f", median(sub$NumberOfIndels,                     na.rm=TRUE)),
            `IQR INDELs`      = sprintf("%.1f", IQR(sub$NumberOfIndels,                        na.rm=TRUE)),
            check.names=FALSE, stringsAsFactors=FALSE
        )
    }))

    n_rows    <- nrow(mean_tbl)
    row_fills <- rep(c("white", "#EEF3FB"), length.out=n_rows)
    tt2 <- ttheme_default(
        core    = list(fg_params=list(cex=0.85),
                       bg_params=list(fill=row_fills, col="gray85")),
        colhead = list(fg_params=list(cex=0.90, fontface="bold", col="white"),
                       bg_params=list(fill="#2166AC", col="#2166AC"))
    )

    grid.newpage()
    pushViewport(viewport(x=0.5, y=0.92, width=1, height=0.12, just=c("centre","centre")))
    grid.text(
        paste0("Per-sample Variant Counts by Filter Stage - Group: ", group_id),
        gp=gpar(fontsize=14, fontface="bold", col="#2166AC")
    )
    popViewport()
    pushViewport(viewport(x=0.5, y=0.83, width=1, height=0.08, just=c("centre","centre")))
    grid.text(
        "Median \u00b1 IQR across all samples",
        gp=gpar(fontsize=9, col="gray40")
    )
    popViewport()
    pushViewport(viewport(x=0.5, y=0.40, width=0.94, height=0.68, just=c("centre","centre")))
    grid.draw(tableGrob(mean_tbl, rows=NULL, theme=tt2))
    popViewport()
    }
}

dev.off()
cat("PDF written:", pdf_file, "\n")
