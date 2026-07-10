suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

thr     <- as.numeric(readLines("threshold.txt")) / 100
tsv_f   <- list.files(".", pattern = "^germline_prevalence_long_", full.names = FALSE)
tsv_f   <- tsv_f[endsWith(tsv_f, ".tsv")]
if (length(tsv_f) == 0L) stop("No germline_prevalence_long_*.tsv found")
df      <- read.table(tsv_f[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Keep only filters that have data
df <- df[!is.na(df[["Prevalence_proportion"]]) & df[["Prevalence_proportion"]] != ".", ]
df[["Prevalence_proportion"]] <- as.numeric(df[["Prevalence_proportion"]])
df <- df[!is.na(df[["Prevalence_proportion"]]), ]

# Ordered factor so facets appear in a logical sequence
filter_order <- c("GERMLINE_FROM_STATS", "VEP_AF_filter", "Bulk_Fail")
df[["Filter"]] <- factor(df[["Filter"]],
                          levels = intersect(filter_order, unique(df[["Filter"]])))

# Pass / Fail relative to the threshold (used for bar fill)
df[["Status"]] <- ifelse(df[["Prevalence_proportion"]] >= thr, "Pass", "Fail")
df[["Status"]] <- factor(df[["Status"]], levels = c("Pass", "Fail"))

# Per-filter summary labels  (n total, n pass, median)
smry <- do.call(rbind, lapply(levels(df[["Filter"]]), function(f) {
    sub   <- df[df[["Filter"]] == f, ]
    n_tot <- nrow(sub)
    n_pas <- sum(sub[["Status"]] == "Pass")
    med   <- median(sub[["Prevalence_proportion"]], na.rm = TRUE)
    data.frame(
        Filter  = f,
        label   = paste0("n = ", format(n_tot, big.mark = ","),
                         "  |  pass = ", format(n_pas, big.mark = ","),
                         "  |  median = ", round(med, 3)),
        stringsAsFactors = FALSE
    )
}))
smry[["Filter"]] <- factor(smry[["Filter"]], levels = levels(df[["Filter"]]))

# Round prevalence to nearest 0.1 so geom_bar produces one visible bar per
# x-axis break (11 bars total: 0.0, 0.1, ..., 1.0).
# Status is determined from the original unrounded value vs the threshold.
df[["Prev_bin"]] <- round(df[["Prevalence_proportion"]], 1L)

p <- ggplot(df, aes(x = Prev_bin, fill = Status)) +
    geom_bar(width = 0.08, colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = thr, linetype = "dashed",
               colour = "black", linewidth = 0.8) +
    geom_text(data = smry,
              aes(x = Inf, y = Inf, label = label),
              hjust = 1.05, vjust = 1.6, size = 3, colour = "grey30",
              inherit.aes = FALSE) +
    facet_wrap(~ Filter, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(Pass = "#A0CC2C", Fail = "#DD14D3"),
                      name = paste0("vs ", round(thr * 100), "% threshold")) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(-0.05, 1.05)) +
    theme_bw(base_size = 11) +
    theme(
        panel.border     = element_rect(colour = "black", linewidth = 0.3),
        strip.background = element_rect(fill = "transparent", colour = NA),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom"
    ) +
    xlab("Proportion of cells with non-REF genotype") +
    ylab("Number of variants") +
    ggtitle(paste0("Germline prevalence distributions — ",
                   gsub("germline_prevalence_long_", "",
                        gsub(".tsv", "", tsv_f[1], fixed = TRUE), fixed = TRUE)))

n_filters <- length(levels(df[["Filter"]]))
ggsave("germline_prevalence_distributions_PLACEHOLDER.png", p,
       width = 9, height = 3.5 * n_filters, dpi = 150, bg = "white")

message("Done.")
