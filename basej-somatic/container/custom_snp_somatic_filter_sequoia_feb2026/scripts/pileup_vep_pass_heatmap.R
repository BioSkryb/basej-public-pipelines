## Heatmap of pileup QC filter status for VEP-PASS variants.
## X-axis : VEP-PASS variants, clustered by failure pattern within two groups.
## Y-axis : single cells (no labels).
## Fill   : which combination of AS / PropClipped / BPPos filters is failing.
## Page 1 : heatmap — variants passing in >= 1 cell
## Page 2 : barplot summary for page 1
## Page 3 : heatmap — variants failing in all cells
## Page 4 : barplot summary for page 3
## Usage  : Rscript pileup_vep_pass_heatmap.R <group_id> <master_table.tsv> <output.pdf>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript pileup_vep_pass_heatmap.R <group_id> <master_table.tsv> <output.pdf>")
group_id    <- args[1]
master_file <- args[2]
output_pdf  <- args[3]

library(ggplot2)
library(scales)
library(gridExtra)

## ── Load master table ─────────────────────────────────────────────────────────
cat("Reading master table...\n")
master <- read.table(master_file, sep = "\t", header = TRUE, quote = "",
                     comment.char = "", stringsAsFactors = FALSE)
vep_pass_ids <- master$VariantId[!is.na(master$VEP_FILTER_STATUS) &
                                  master$VEP_FILTER_STATUS == "PASS"]
cat("VEP PASS variants:", length(vep_pass_ids), "\n")

## ── Load pileup files ─────────────────────────────────────────────────────────
pileup_files <- list.files(".", pattern = paste0("^res_pileup_all_group_", group_id,
                                                   "_chr.*[.]tsv$"), full.names = TRUE)
if (length(pileup_files) == 0) stop("No res_pileup files found for group ", group_id)
cat("Loading", length(pileup_files), "pileup files...\n")
pil <- do.call(rbind, lapply(pileup_files, read.table, sep = "\t", header = TRUE,
                              quote = "", comment.char = "", stringsAsFactors = FALSE))

pil <- pil[pil$ALT != "REF", ]
pil <- pil[pil$VariantId %in% vep_pass_ids, ]
cat("VEP PASS variants found in pileup:", length(unique(pil$VariantId)), "\n")
cat("Samples:", length(unique(pil$SampleId)), "\n")
cat("Variant x Sample pairs:", nrow(pil), "\n\n")

## ── Filter status per row ─────────────────────────────────────────────────────
pil$as_fail   <- pil$AS_Filter          == "Fail"
pil$clip_fail <- pil$PropClipped_Filter == "Fail"
pil$bp_fail   <- pil$BPPos_Filter       == "Fail"

pil$filter_status <- "Pass"
pil$filter_status[  pil$as_fail & !pil$clip_fail & !pil$bp_fail] <- "AS"
pil$filter_status[ !pil$as_fail &  pil$clip_fail & !pil$bp_fail] <- "Clipping"
pil$filter_status[ !pil$as_fail & !pil$clip_fail &  pil$bp_fail] <- "BPPos"
pil$filter_status[  pil$as_fail &  pil$clip_fail & !pil$bp_fail] <- "AS + Clipping"
pil$filter_status[  pil$as_fail & !pil$clip_fail &  pil$bp_fail] <- "AS + BPPos"
pil$filter_status[ !pil$as_fail &  pil$clip_fail &  pil$bp_fail] <- "Clipping + BPPos"
pil$filter_status[  pil$as_fail &  pil$clip_fail &  pil$bp_fail] <- "All three fail"

cat("Filter status distribution:\n")
print(sort(table(pil$filter_status), decreasing = TRUE)); cat("\n")

## ── Split variants ────────────────────────────────────────────────────────────
pass_per_var <- tapply(pil$filter_status == "Pass", pil$VariantId, sum, na.rm = TRUE)
has_any_pass <- names(pass_per_var)[pass_per_var >= 1]
never_pass   <- names(pass_per_var)[pass_per_var == 0]
cat("Variants with >= 1 Pass cell:", length(has_any_pass), "\n")
cat("Variants with 0  Pass cells: ", length(never_pass),   "\n\n")

## ── Clustering ────────────────────────────────────────────────────────────────
status_num <- c("Pass" = 0, "No data" = 0,
                "AS" = 1, "Clipping" = 2, "BPPos" = 3,
                "AS + Clipping" = 4, "AS + BPPos" = 5,
                "Clipping + BPPos" = 6, "All three fail" = 7)

all_samples  <- sort(unique(pil$SampleId))
all_variants <- unique(pil$VariantId)

mat <- matrix(0, nrow = length(all_variants), ncol = length(all_samples),
              dimnames = list(all_variants, all_samples))
for (i in seq_len(nrow(pil)))
    mat[pil$VariantId[i], pil$SampleId[i]] <- status_num[pil$filter_status[i]]

cluster_order <- function(ids) {
    if (length(ids) <= 1) return(ids)
    hc <- hclust(dist(mat[ids, , drop = FALSE], method = "euclidean"), method = "ward.D2")
    ids[hc$order]
}

ordered_has_pass <- cluster_order(has_any_pass)
ordered_never    <- cluster_order(never_pass)

## ── Colour palette ────────────────────────────────────────────────────────────
status_levels <- c("Pass", "AS", "Clipping", "BPPos",
                   "AS + Clipping", "AS + BPPos", "Clipping + BPPos",
                   "All three fail", "No data")
pal <- c(
    "Pass"             = "#CCCCCC",
    "AS"               = "#E6AB02",
    "Clipping"         = "#377EB8",
    "BPPos"            = "#E41A1C",
    "AS + Clipping"    = "#A65628",
    "AS + BPPos"       = "#984EA3",
    "Clipping + BPPos" = "#FF7F00",
    "All three fail"   = "#111111",
    "No data"          = "#F5F5F5"
)

## ── Build plot_df ─────────────────────────────────────────────────────────────
make_plot_df <- function(var_ids) {
    fg <- expand.grid(VariantId = var_ids, SampleId = all_samples,
                      stringsAsFactors = FALSE)
    df <- merge(fg, pil[, c("VariantId","SampleId","filter_status")],
                by = c("VariantId","SampleId"), all.x = TRUE)
    df$filter_status[is.na(df$filter_status)] <- "No data"
    df$filter_status <- factor(df$filter_status, levels = status_levels)
    df$VariantId     <- factor(df$VariantId, levels = var_ids)
    df$SampleId      <- factor(df$SampleId,  levels = all_samples)
    df$var_idx       <- as.integer(df$VariantId)
    df$smp_idx       <- as.integer(df$SampleId)
    df
}

df1 <- make_plot_df(ordered_has_pass)
df2 <- make_plot_df(ordered_never)

## ── Per-group statistics ──────────────────────────────────────────────────────
group_stats <- function(var_ids, df_plot) {
    d            <- df_plot[df_plot$filter_status != "No data", ]
    n_pairs_data <- nrow(d)
    n_pass       <- sum(d$filter_status == "Pass")
    n_fail       <- n_pairs_data - n_pass
    list(n_vars       = length(var_ids),
         n_pairs_data = n_pairs_data,
         n_pass       = n_pass,
         n_fail       = n_fail)
}
s1 <- group_stats(ordered_has_pass, df1)
s2 <- group_stats(ordered_never,    df2)

## ── Shared theme ──────────────────────────────────────────────────────────────
theme_heatmap <- theme_bw(base_size = 9) +
    theme(
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_text(size = 8, face = "bold"),
        legend.text      = element_text(size = 8),
        plot.title       = element_text(face = "bold", size = 10),
        plot.subtitle    = element_text(size = 8)
    )

legend_guide <- scale_fill_manual(
    values = pal, name = "Filter failing", drop = FALSE,
    guide  = guide_legend(nrow = 2, byrow = TRUE,
                          override.aes = list(size = 5))
)

make_heatmap <- function(df, title, subtitle) {
    ggplot(df, aes(x = var_idx, y = smp_idx, fill = filter_status)) +
        geom_tile() +
        legend_guide +
        labs(title = title, subtitle = subtitle,
             x = paste0("Variants (n=", nlevels(df$VariantId),
                        ")  —  clustered by failure pattern"),
             y = paste0("Cells (n=", length(all_samples), ")")) +
        theme_heatmap
}

## ── Barplot: proportion of pairs per filter status ───────────────────────────
make_barplot <- function(df, title, stats) {
    # Exclude "No data" from proportion calculation
    d <- df[df$filter_status != "No data", ]
    total <- nrow(d)
    cnt   <- as.data.frame(table(filter_status = d$filter_status))
    cnt   <- cnt[cnt$Freq > 0, ]
    cnt$prop  <- cnt$Freq / total
    cnt$label <- paste0(cnt$Freq, "\n(", round(100 * cnt$prop, 1), "%)")
    cnt$filter_status <- factor(cnt$filter_status, levels = status_levels)

    sub <- paste0(
        stats$n_vars, " unique VEP-PASS variants in pileup  |  ",
        stats$n_pairs_data, " variant\u00d7cell pairs with data  ",
        "(Pass: ", stats$n_pass, "  [", round(100 * stats$n_pass / stats$n_pairs_data, 1), "%],  ",
        "Fail: ", stats$n_fail, "  [", round(100 * stats$n_fail / stats$n_pairs_data, 1), "%])"
    )

    ggplot(cnt, aes(x = filter_status, y = prop, fill = filter_status)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = label), vjust = -0.3, size = 3) +
        scale_fill_manual(values = pal, drop = FALSE, guide = "none") +
        scale_y_continuous(labels = percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.15))) +
        labs(title    = title,
             subtitle = sub,
             x = NULL, y = "Proportion of pairs") +
        theme_bw(base_size = 10) +
        theme(axis.text.x  = element_text(angle = 35, hjust = 1, size = 9),
              panel.grid.major.x = element_blank(),
              plot.title   = element_text(face = "bold", size = 10),
              plot.subtitle = element_text(size = 8))
}

p1 <- make_heatmap(df1,
    title    = paste0("Heatmap — VEP-PASS variants with pileup data, passing in \u2265 1 cell  |  Group: ", group_id),
    subtitle = paste0(
        s1$n_vars, " unique VEP-PASS variants in pileup  |  ",
        s1$n_pairs_data, " variant\u00d7cell pairs with data  ",
        "(Pass: ", s1$n_pass, "  [", round(100 * s1$n_pass / s1$n_pairs_data, 1), "%],  ",
        "Fail: ", s1$n_fail, "  [", round(100 * s1$n_fail / s1$n_pairs_data, 1), "%])\n",
        "Clustered by failure profile (ward.D2).  Grey = Pass.  Near-white = No coverage."
    ))

b1 <- make_barplot(df1,
    title = paste0("Filter status proportions — VEP-PASS variants passing in \u2265 1 cell  |  Group: ", group_id),
    stats = s1)

p2 <- make_heatmap(df2,
    title    = paste0("Heatmap — VEP-PASS variants with pileup data, failing in ALL cells  |  Group: ", group_id),
    subtitle = paste0(
        s2$n_vars, " unique VEP-PASS variants in pileup  |  ",
        s2$n_pairs_data, " variant\u00d7cell pairs with data  ",
        "(Pass: ", s2$n_pass, "  [", round(100 * s2$n_pass / s2$n_pairs_data, 1), "%],  ",
        "Fail: ", s2$n_fail, "  [", round(100 * s2$n_fail / s2$n_pairs_data, 1), "%])\n",
        "Clustered by failure profile (ward.D2).  These variants never pass pileup QC in any cell."
    ))

b2 <- make_barplot(df2,
    title = paste0("Filter status proportions — VEP-PASS variants failing in all cells  |  Group: ", group_id),
    stats = s2)

## ── Write PDF ─────────────────────────────────────────────────────────────────
pdf(output_pdf, width = 22, height = 8)
print(p1)
print(b1)
print(p2)
print(b2)
dev.off()
cat("Written:", output_pdf, "(4 pages)\n")
