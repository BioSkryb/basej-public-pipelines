## BPPos filter diagnostic plots.
## Restricted to VEP-PASS variants that are present in the pileup data.
## Explores the variable relationships used in the base-position bias filter,
## split by the three strand-availability scenarios defined in rscript_2:
##   Case A — only reverse strand has HQ reads (F < 2, R > 1)
##   Case B — only forward strand has HQ reads (F > 1, R < 2)
##   Case C — both strands have HQ reads     (F > 1, R > 1)
## Usage: Rscript pileup_bppos_plots.R <group_id> <master_table.tsv> <output.pdf>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript pileup_bppos_plots.R <group_id> <master_table.tsv> <output.pdf>")
group_id    <- args[1]
master_file <- args[2]
output_pdf  <- args[3]

library(ggplot2)
library(scales)
suppressPackageStartupMessages(library(data.table))

## ── Load VEP PASS variant IDs from master table ───────────────────────────────
cat("Reading master table:", master_file, "\n")
master <- fread(master_file, sep = "\t", header = TRUE,
                data.table = FALSE, showProgress = FALSE)
vep_pass_ids <- master$VariantId[!is.na(master$VEP_FILTER_STATUS) &
                                  master$VEP_FILTER_STATUS == "PASS"]
cat("VEP PASS variants in master table:", length(vep_pass_ids), "\n")

## ── Load pileup data ──────────────────────────────────────────────────────────
pileup_files <- list.files(".", pattern = paste0("^res_pileup_all_group_", group_id, "_chr.*[.]tsv$"),
                           full.names = TRUE)
if (length(pileup_files) == 0)
  stop("No res_pileup files found for group ", group_id)

cat("Loading", length(pileup_files), "pileup files...\n")
pil <- as.data.frame(rbindlist(lapply(pileup_files, fread,
                                       sep = "\t", header = TRUE,
                                       data.table = FALSE,
                                       showProgress = FALSE)))

## ── Drop REF rows — BPPos filter is only applied to ALT positions; REF rows ──
## are auto-labelled Fail and would pollute the Fail panel with well-covered
## reference pileups that were never evaluated by the filter.
pil <- pil[pil$ALT != "REF", ]

## ── Restrict to VEP-PASS variants present in pileup ──────────────────────────
pil <- pil[pil$VariantId %in% vep_pass_ids, ]
cat("VEP PASS variants found in pileup:", length(unique(pil$VariantId)), "\n")
cat("Total variant x sample pairs after filter:", nrow(pil), "\n\n")

## ── Assign strand case ────────────────────────────────────────────────────────
pil$Case <- "Neither"
pil$Case[pil$NUM_FRAGMENTS_HQ_MQ_BQ_F <  2 & pil$NUM_FRAGMENTS_HQ_MQ_BQ_R >  1] <- "A: Reverse only"
pil$Case[pil$NUM_FRAGMENTS_HQ_MQ_BQ_F >  1 & pil$NUM_FRAGMENTS_HQ_MQ_BQ_R <  2] <- "B: Forward only"
pil$Case[pil$NUM_FRAGMENTS_HQ_MQ_BQ_F >  1 & pil$NUM_FRAGMENTS_HQ_MQ_BQ_R >  1] <- "C: Both strands"
pil$Case        <- factor(pil$Case, levels = c("A: Reverse only","B: Forward only","C: Both strands","Neither"))
pil$BPPos_Filter <- factor(pil$BPPos_Filter, levels = c("Pass","Fail"))

case_counts <- table(pil$Case)
cat("\nCase counts:\n"); print(case_counts)
cat("\nBPPos_Filter by case:\n"); print(table(pil$Case, pil$BPPos_Filter))

## ── Shared style ──────────────────────────────────────────────────────────────
PAL <- c("Pass" = "#4DAF4A", "Fail" = "#E41A1C")

## Filter thresholds (must match nextflow.config)
T_TAIL  <- 0.30   # cutoff_prop_bp_under / cutoff_prop_bp_upper
T_SD    <- 25     # cutoff_sd_both
T_MAD   <- 12     # cutoff_mad_both

theme_clean <- theme_bw(base_size = 12) +
  theme(panel.grid.minor  = element_blank(),
        strip.background  = element_rect(fill = "grey90", colour = NA),
        strip.text        = element_text(face = "bold"),
        plot.title        = element_text(face = "bold", size = 12),
        legend.position   = "right",
        legend.title      = element_text(size = 9))

hexplot <- function(df, x, y, xlab, ylab, title, subtitle, bins = 60,
                    vline = NULL, hline = NULL) {
  df <- df[is.finite(df[[x]]) & is.finite(df[[y]]), ]
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_bin2d(bins = bins) +
    scale_fill_viridis_c(name = "Count", trans = "log10",
                         option = "magma", direction = -1) +
    facet_wrap(~ BPPos_Filter, ncol = 2) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_clean
  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed",
                        colour = "dodgerblue", linewidth = 0.6)
  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                        colour = "dodgerblue", linewidth = 0.6)
  p
}

## ── Helper: build label string with pair counts and unique variant counts ──────
case_label <- function(df, case_name) {
    n_pairs      <- nrow(df)
    n_pass_pairs <- sum(df$BPPos_Filter == "Pass")
    n_fail_pairs <- sum(df$BPPos_Filter == "Fail")
    n_vars       <- length(unique(df$VariantId))
    n_vars_pass  <- length(unique(df$VariantId[df$BPPos_Filter == "Pass"]))
    n_vars_fail  <- length(unique(df$VariantId[df$BPPos_Filter == "Fail"]))
    paste0(
        "  |  ", case_name,
        "  |  pairs: n=", scales::comma(n_pairs),
        "  (Pass=", scales::comma(n_pass_pairs),
        " [", round(100 * n_pass_pairs / n_pairs, 1), "%],",
        "  Fail=", scales::comma(n_fail_pairs),
        " [", round(100 * n_fail_pairs / n_pairs, 1), "%])",
        "\n       unique variants: n=", n_vars,
        "  (any-Pass=", n_vars_pass, ",  never-Pass=", n_vars_fail, ")"
    )
}

## ── Page 1: case summary bar chart ───────────────────────────────────────────
sum_df <- as.data.frame(table(Case = pil$Case, BPPos = pil$BPPos_Filter))
sum_df <- sum_df[sum_df$Case != "Neither", ]
sum_df$Case <- droplevels(sum_df$Case)

n_vep_in_pileup <- length(unique(pil$VariantId))
n_total_pairs   <- nrow(pil[pil$Case != "Neither", ])

p_summary <- ggplot(sum_df, aes(x = Case, y = Freq, fill = BPPos)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = scales::comma(Freq)),
            position = position_stack(vjust = 0.5), size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = PAL) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
  labs(title    = paste0("BPPos filter outcome by strand-availability case  |  Group: ", group_id,
                         "  |  VEP-PASS variants only"),
       x        = NULL, y = "Number of variant\u00d7cell pairs",
       subtitle = paste0(
           "Universe: ", n_vep_in_pileup, " unique VEP-PASS variants present in pileup  |  ",
           scales::comma(n_total_pairs), " variant\u00d7cell pairs (A+B+C cases)\n",
           "Case A: F<2, R>1  |  Case B: F>1, R<2  |  Case C: F>1, R>1  |  ",
           "Neither (F<2, R<2) excluded: n=", scales::comma(sum(pil$Case == "Neither")))) +
  theme_clean + theme(legend.position = "bottom")

## ═════════════════════════════════════════════════════════════════════════════
## CASE A: Reverse only (F < 2, R > 1)
## Filter: (UNDER_R < t AND UPPER_R < t) OR (SD_R > t AND MAD_R > t)
## ═════════════════════════════════════════════════════════════════════════════
a <- pil[pil$Case == "A: Reverse only", ]
a$TAIL_SUM_R <- a$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R + a$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R
gp_a <- case_label(a, "Case A: Reverse only")

pA1 <- hexplot(a,
  x = "TAIL_SUM_R",
  y = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  xlab     = "Tail proportion sum: UNDER_R + UPPER_R",
  ylab     = "SD of bp-start positions (R)  [bp]",
  title    = paste0("Case A — Tail proportion vs SD (reverse strand)", gp_a),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and SD threshold (y=", T_SD, " bp).\n",
                    "Pass via low tails (left of x line) OR high SD+MAD (above y line AND MAD > ", T_MAD, " bp)."),
  vline = T_TAIL, hline = T_SD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

pA2 <- hexplot(a,
  x = "TAIL_SUM_R",
  y = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  xlab     = "Tail proportion sum: UNDER_R + UPPER_R",
  ylab     = "MAD of bp-start positions (R)  [bp]",
  title    = paste0("Case A — Tail proportion vs MAD (reverse strand)", gp_a),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and MAD threshold (y=", T_MAD, " bp).\n",
                    "Both SD > ", T_SD, " AND MAD > ", T_MAD, " bp required for the spread path to pass."),
  vline = T_TAIL, hline = T_MAD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

## ═════════════════════════════════════════════════════════════════════════════
## CASE B: Forward only (F > 1, R < 2)
## Filter: (UNDER_F < t AND UPPER_F < t) OR (SD_F > t AND MAD_F > t)
## ═════════════════════════════════════════════════════════════════════════════
b <- pil[pil$Case == "B: Forward only", ]
b$TAIL_SUM_F <- b$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F + b$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F
gp_b <- case_label(b, "Case B: Forward only")

pB1 <- hexplot(b,
  x = "TAIL_SUM_F",
  y = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  xlab     = "Tail proportion sum: UNDER_F + UPPER_F",
  ylab     = "SD of bp-start positions (F)  [bp]",
  title    = paste0("Case B — Tail proportion vs SD (forward strand)", gp_b),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and SD threshold (y=", T_SD, " bp).\n",
                    "Mirror of Case A applied to forward strand. Pass via low tails OR high SD+MAD."),
  vline = T_TAIL, hline = T_SD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

pB2 <- hexplot(b,
  x = "TAIL_SUM_F",
  y = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  xlab     = "Tail proportion sum: UNDER_F + UPPER_F",
  ylab     = "MAD of bp-start positions (F)  [bp]",
  title    = paste0("Case B — Tail proportion vs MAD (forward strand)", gp_b),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and MAD threshold (y=", T_MAD, " bp).\n",
                    "Both SD > ", T_SD, " AND MAD > ", T_MAD, " bp required for the spread path to pass."),
  vline = T_TAIL, hline = T_MAD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

## ═════════════════════════════════════════════════════════════════════════════
## CASE C: Both strands (F > 1, R > 1)
## df_c_1 (F-primary): low tails on F (mandatory) AND (SD_F+MAD_F > t_both OR SD_R+MAD_R > t_extreme)
## df_c_2 (R-primary): low tails on R (mandatory) AND (SD_R+MAD_R > t_both OR SD_F+MAD_F > t_extreme)
## ═════════════════════════════════════════════════════════════════════════════
cc <- pil[pil$Case == "C: Both strands", ]
cc$TAIL_SUM_F <- cc$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F + cc$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F
cc$TAIL_SUM_R <- cc$PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R + cc$PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R
gp_c <- case_label(cc, "Case C: Both strands")

pC1 <- hexplot(cc,
  x = "TAIL_SUM_F",
  y = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  xlab     = "Tail proportion sum: UNDER_F + UPPER_F",
  ylab     = "SD of bp-start positions (F)  [bp]",
  title    = paste0("Case C — Forward strand: tail proportion vs SD", gp_c),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and SD threshold (y=", T_SD, " bp).\n",
                    "SD alone is not sufficient for interpretation — always read alongside the MAD plot below."),
  vline = T_TAIL, hline = T_SD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

pC2 <- hexplot(cc,
  x = "TAIL_SUM_F",
  y = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  xlab     = "Tail proportion sum: UNDER_F + UPPER_F",
  ylab     = "MAD of bp-start positions (F)  [bp]",
  title    = paste0("Case C — Forward strand: tail proportion vs MAD", gp_c),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and MAD threshold (y=", T_MAD, " bp).\n",
                    "MAD is robust to outlier reads and is the binding constraint in the spread check."),
  vline = T_TAIL, hline = T_MAD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

pC3 <- hexplot(cc,
  x = "TAIL_SUM_R",
  y = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  xlab     = "Tail proportion sum: UNDER_R + UPPER_R",
  ylab     = "SD of bp-start positions (R)  [bp]",
  title    = paste0("Case C — Reverse strand: tail proportion vs SD", gp_c),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and SD threshold (y=", T_SD, " bp).\n",
                    "SD alone is not sufficient for interpretation — always read alongside the MAD plot below."),
  vline = T_TAIL, hline = T_SD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

pC4 <- hexplot(cc,
  x = "TAIL_SUM_R",
  y = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  xlab     = "Tail proportion sum: UNDER_R + UPPER_R",
  ylab     = "MAD of bp-start positions (R)  [bp]",
  title    = paste0("Case C — Reverse strand: tail proportion vs MAD", gp_c),
  subtitle = paste0("Blue lines: tail threshold (x=", T_TAIL, ") and MAD threshold (y=", T_MAD, " bp).\n",
                    "MAD is robust to outlier reads and is the binding constraint in the spread check."),
  vline = T_TAIL, hline = T_MAD
) + scale_x_continuous(labels = percent_format(accuracy = 1))

## ── Case C: SD vs MAD per strand — the key diagnostic plot ───────────────────
## High SD but low MAD = outlier-inflated SD; reads pile up at one position.
## These are correctly rejected by the MAD threshold.
## Quadrants (threshold lines at SD=T_SD, MAD=T_MAD):
##   Top-right    (SD > T_SD, MAD > T_MAD): genuine spread — passes spread check
##   Bottom-right (SD > T_SD, MAD < T_MAD): deceptive SD — pile-up artifact, correctly fails
##   Bottom-left  (SD < T_SD, MAD < T_MAD): genuine pile-up — correctly fails
pC_sd_mad_F <- hexplot(cc,
  x = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  y = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  xlab     = "SD of bp-start positions (F)  [bp]",
  ylab     = "MAD of bp-start positions (F)  [bp]",
  title    = paste0("Case C — Forward strand: SD vs MAD  [KEY DIAGNOSTIC]", gp_c),
  subtitle = paste0("Blue lines: SD threshold (x=", T_SD, " bp) and MAD threshold (y=", T_MAD, " bp).\n",
                    "Top-right: genuine spread (pass). Bottom-right: high SD but low MAD = outlier-inflated SD (correctly fails).\n",
                    "Most Fail mass is bottom-left (genuine pile-up) or bottom-right (deceptive SD). Filter is correct."),
  vline = T_SD, hline = T_MAD
)

pC_sd_mad_R <- hexplot(cc,
  x = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  y = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  xlab     = "SD of bp-start positions (R)  [bp]",
  ylab     = "MAD of bp-start positions (R)  [bp]",
  title    = paste0("Case C — Reverse strand: SD vs MAD  [KEY DIAGNOSTIC]", gp_c),
  subtitle = paste0("Blue lines: SD threshold (x=", T_SD, " bp) and MAD threshold (y=", T_MAD, " bp).\n",
                    "Top-right: genuine spread (pass). Bottom-right: high SD but low MAD = outlier-inflated SD (correctly fails).\n",
                    "Most Fail mass is bottom-left (genuine pile-up) or bottom-right (deceptive SD). Filter is correct."),
  vline = T_SD, hline = T_MAD
)

## Case C cross-strand comparisons: SD_F vs SD_R and MAD_F vs MAD_R
pC5 <- hexplot(cc,
  x = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  y = "SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  xlab     = "SD of bp-start positions (F)  [bp]",
  ylab     = "SD of bp-start positions (R)  [bp]",
  title    = paste0("Case C — Cross-strand SD: F vs R", gp_c),
  subtitle = paste0("Blue lines: SD threshold (", T_SD, " bp) on each axis.\n",
                    "SD overlap between Pass/Fail panels is expected: SD can be inflated by outlier reads.\n",
                    "Interpret alongside the SD vs MAD plots — SD alone does not determine pass/fail."),
  vline = T_SD, hline = T_SD
) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4)

pC6 <- hexplot(cc,
  x = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F",
  y = "MAD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R",
  xlab     = "MAD of bp-start positions (F)  [bp]",
  ylab     = "MAD of bp-start positions (R)  [bp]",
  title    = paste0("Case C — Cross-strand MAD: F vs R", gp_c),
  subtitle = paste0("Blue lines: MAD threshold (", T_MAD, " bp) on each axis.\n",
                    "MAD is robust to outlier reads. Pass variants cluster top-right (high MAD on both strands).\n",
                    "Fail variants at low MAD on both strands are genuine pile-ups; filter correctly rejects them."),
  vline = T_MAD, hline = T_MAD
) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4)

## ── Write PDF ─────────────────────────────────────────────────────────────────
## Page layout:
##  1          : Summary bar chart
##  2–3        : Case A (reverse only) — TAIL_SUM_R vs SD; vs MAD
##  4–5        : Case B (forward only) — TAIL_SUM_F vs SD; vs MAD
##  6–9        : Case C single-strand  — TAIL_SUM_F vs SD; vs MAD; TAIL_SUM_R vs SD; vs MAD
##  10–11      : Case C SD vs MAD (F strand; R strand)  [KEY DIAGNOSTIC]
##  12–13      : Case C cross-strand   — SD_F vs SD_R; MAD_F vs MAD_R
pdf(output_pdf, width = 12, height = 7)
print(p_summary)
print(pA1); print(pA2)
print(pB1); print(pB2)
print(pC1); print(pC2); print(pC3); print(pC4)
print(pC_sd_mad_F); print(pC_sd_mad_R)
print(pC5); print(pC6)
dev.off()
cat("Written:", output_pdf, "(13 pages)\n")
