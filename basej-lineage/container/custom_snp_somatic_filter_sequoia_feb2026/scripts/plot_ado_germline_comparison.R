suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(grid)
    library(RColorBrewer)
})

# ── Interval metadata ────────────────────────────────────────────────────────
# Levels as stored in merged_ADO (brackets stripped by CONCAT_SUMMARY_ADO_INTERVALS)
interval_levels <- c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",
                     "0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1")
# Display labels matching the original ADO comparison figure style
interval_labels <- c("[0-0.1)","[0.1-0.2)","[0.2-0.3)","[0.3-0.4)","[0.4-0.5)",
                     "[0.5-0.6)","[0.6-0.7)","[0.7-0.8)","[0.8-0.9)","[0.9-1]")
prov_order <- c("stats", "vep", "bulk")

# ── Helper: extract provenance from filename without regex back-refs ──────────
# "merged_ADO_stats.tsv"         -> "stats"   (prefix len = 11, suffix len = 4)
# "merged_ADO_summary_stats.tsv" -> "stats"   (prefix len = 19, suffix len = 4)
prov_from_ado     <- function(f) substr(basename(f), 12L, nchar(basename(f)) - 4L)
prov_from_summary <- function(f) substr(basename(f), 20L, nchar(basename(f)) - 4L)

# ── 1. Load merged_ADO_*.tsv ─────────────────────────────────────────────────
all_f     <- list.files(".", full.names = FALSE)
ado_files <- all_f[startsWith(all_f, "merged_ADO_") &
                   endsWith(all_f, ".tsv") &
                   !grepl("summary", all_f, fixed = TRUE)]
if (length(ado_files) == 0L) stop("No merged_ADO_*.tsv files found")

df_ado <- do.call(rbind, lapply(ado_files, function(f) {
    prov <- prov_from_ado(f)
    dat  <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # File_Interval format: "{SampleId}_{provenance}_{interval}"
    # Use fixed-string search so no regex metacharacter issues
    sep <- paste0("_", prov, "_")
    pos <- regexpr(sep, dat[["File_Interval"]], fixed = TRUE)
    dat[["SampleId"]]   <- substr(dat[["File_Interval"]], 1L, pos - 1L)
    dat[["Interval"]]   <- substr(dat[["File_Interval"]],
                                  pos + nchar(sep),
                                  nchar(dat[["File_Interval"]]))
    dat[["Provenance"]] <- prov
    dat
}))

df_ado[["Interval"]]   <- factor(df_ado[["Interval"]],
                                  levels = interval_levels,
                                  labels = interval_labels)
df_ado[["Provenance"]] <- factor(df_ado[["Provenance"]],
                                  levels = intersect(prov_order,
                                                     unique(df_ado[["Provenance"]])))

# ── 2. Load merged_ADO_summary_*.tsv ─────────────────────────────────────────
sum_files <- all_f[startsWith(all_f, "merged_ADO_summary_") & endsWith(all_f, ".tsv")]
if (length(sum_files) == 0L) stop("No merged_ADO_summary_*.tsv files found")

df_sum <- do.call(rbind, lapply(sum_files, function(f) {
    prov <- prov_from_summary(f)
    dat  <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    dat[["Provenance"]] <- prov
    # sample_name: "{SampleId}_{provenance}" — strip trailing "_{provenance}"
    sn <- dat[["sample_name"]]
    dat[["SampleId"]] <- substr(sn, 1L, nchar(sn) - nchar(prov) - 1L)
    dat
}))

df_sum[["Provenance"]] <- factor(df_sum[["Provenance"]],
                                  levels = intersect(prov_order,
                                                     unique(df_sum[["Provenance"]])))

# ── Scalable colour palette for connecting lines in the summary panel ─────────
all_samples <- union(unique(df_ado[["SampleId"]]), unique(df_sum[["SampleId"]]))
n_samp      <- length(all_samples)
if (n_samp <= 8L) {
    pal <- brewer.pal(max(3L, n_samp), "Set1")[seq_len(n_samp)]
} else {
    pal <- colorRampPalette(brewer.pal(8L, "Set1"))(n_samp)
}
names(pal) <- all_samples

# ── Per-provenance annotation: "sites × samples = N total" ───────────────────
# For each provenance compute n_sites_per_sample, n_samples, and grand total.
make_count_label <- function(prov_factor, df) {
    do.call(rbind, lapply(levels(prov_factor), function(p) {
        sub     <- df[as.character(df[["Provenance"]]) == p, ]
        n_s     <- length(unique(sub[["SampleId"]]))
        n_tot   <- sum(sub[["Freq"]])
        n_sites <- if (n_s > 0L) round(n_tot / n_s) else 0L
        data.frame(
            Provenance = p,
            label      = paste0(format(n_sites, big.mark = ",", trim = TRUE),
                                " \u00d7 ", n_s, " = ",
                                format(n_tot, big.mark = ",", trim = TRUE)),
            stringsAsFactors = FALSE
        )
    }))
}
ado_vcounts <- make_count_label(df_ado[["Provenance"]], df_ado)
ado_vcounts[["Provenance"]] <- factor(ado_vcounts[["Provenance"]],
                                       levels = levels(df_ado[["Provenance"]]))
# x placement: last discrete level (Inf not valid on discrete x scales in ggplot2 3.4)
ado_vcounts[["x_lab"]] <- factor(interval_labels[length(interval_labels)],
                                  levels = interval_labels)

# Same label for the summary panel (sites from df_ado, samples from df_sum)
sum_vcounts <- make_count_label(df_ado[["Provenance"]], df_ado)
sum_vcounts[["Provenance"]] <- factor(sum_vcounts[["Provenance"]],
                                       levels = levels(df_sum[["Provenance"]]))

# ── 3. Distribution plot (left panel) ────────────────────────────────────────
# Shading for the balanced-allele zone [0.2-0.8]: interval labels 3-8.
# annotate("rect") conflicts with discrete x in ggplot2 3.4.x — use geom_tile.
shade_df <- data.frame(
    Interval = factor(interval_labels[3:8], levels = interval_labels),
    y        = 0.25
)
p_dist <- ggplot(df_ado, aes(x = Interval, y = Prop)) +
    geom_tile(data = shade_df, aes(x = Interval, y = y, width = 1, height = Inf),
              fill = "#D9D9D9", alpha = 0.5, inherit.aes = FALSE) +
    geom_point(size = 1.8, alpha = 0.6, colour = "grey40") +
    geom_line(aes(group = SampleId, colour = SampleId), alpha = 0.55, linewidth = 0.4) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = 1), colour = "black", linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point",
                 aes(group = 1), shape = 15, colour = "red", size = 3.5) +
    geom_text(data = ado_vcounts,
              aes(x = x_lab, y = Inf, label = label),
              hjust = 1.05, vjust = 1.5, size = 3, colour = "grey30",
              inherit.aes = FALSE) +
    facet_grid(. ~ Provenance, scales = "free_y") +
    scale_colour_manual(values = pal) +
    theme_bw(base_size = 10) +
    theme(
        panel.border     = element_rect(colour = "black", linewidth = 0.3),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position  = "none",
        strip.background = element_rect(fill = "transparent", colour = NA)
    ) +
    ylab("Proportion of total ground truth sites") +
    xlab(NULL) +
    scale_x_discrete(limits = interval_labels) +
    scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0, 1)) +
    ggtitle("Allelic balance: Distributions  (per germline filter set)")

# ── 4. Summary / ADO-proportion plot (right panel) ───────────────────────────
# Boxplot + uncoloured points; coloured dashed lines connect the same cell.
p_sum <- ggplot(df_sum, aes(x = Provenance, y = ADO_PERC)) +
    geom_hline(yintercept = 0.8, colour = "red", linewidth = 0.7) +
    geom_boxplot(outlier.shape = NA, width = 0.35, colour = "black", fill = NA) +
    geom_line(aes(group = SampleId, colour = SampleId),
              linetype = "dashed", alpha = 0.55) +
    geom_point(size = 3, colour = "grey40") +
    geom_text(data = sum_vcounts,
              aes(x = Provenance, y = Inf, label = label),
              vjust = 1.5, size = 2.5, colour = "grey30",
              inherit.aes = FALSE) +
    scale_colour_manual(values = pal) +
    theme_bw(base_size = 10) +
    theme(
        panel.border     = element_rect(colour = "black", linewidth = 0.3),
        legend.position  = "none",
        strip.background = element_rect(fill = "grey92")
    ) +
    ylab("Proportion in allelic balance [0.2-0.8]") +
    xlab("Germline filter set") +
    scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0, 1)) +
    ggtitle("Allelic balance: Summary\n(lines = same cell)")

# ── 5. Save individual plots ──────────────────────────────────────────────────
ggsave("ADO_germline_dist.png",    p_dist,
       width = 13, height = 6, dpi = 150, bg = "white")
ggsave("ADO_germline_summary.png", p_sum,
       width = 5,  height = 6, dpi = 150, bg = "white")

# ── 6. Combine side-by-side with grid viewports ───────────────────────────────
# Distribution (72 %) left, summary (28 %) right — mirrors the original figure.
png("ADO_germline_comparison.png",
    width = 2700L, height = 900L, res = 150L, bg = "white")
grid.newpage()
pushViewport(viewport(x = 0, y = 0, width = 0.72, height = 1,
                      just = c("left", "bottom")))
print(p_dist, newpage = FALSE)
popViewport()
pushViewport(viewport(x = 0.72, y = 0, width = 0.28, height = 1,
                      just = c("left", "bottom")))
print(p_sum, newpage = FALSE)
popViewport()
dev.off()

write.table(df_sum, "ADO_germline_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Done: ADO_germline_dist.png  ADO_germline_summary.png  ADO_germline_comparison.png  ADO_germline_summary.tsv")
