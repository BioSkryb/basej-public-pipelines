## Binomial / Betabinomial Filter Plots (ggplot2).
## Called from CUSTOM_VARIANT_FILTER_PROVENANCE.
## Usage: Rscript binom_ggplots.R <master_table.tsv> <group_id> <output.pdf>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript binom_ggplots.R <master_table.tsv> <group_id> <output.pdf>")
master_tsv  <- args[1]
group_id    <- args[2]
output_pdf  <- args[3]

library(ggplot2)
library(paletteer)
library(scales)

d <- read.table(master_tsv, sep = "\t", header = TRUE, quote = "", comment.char = "")

## Verdict labels & colours
verdict_map <- c(
  "TRUE_TRUE_TRUE"   = "Somatic (Pass)",
  "TRUE_FALSE_FALSE" = "Germline: both",
  "TRUE_FALSE_TRUE"  = "Germline: Binom only",
  "TRUE_TRUE_FALSE"  = "Germline: Betabinom only"
)
pal <- setNames(as.character(paletteer_d("RColorBrewer::Set1")[1:4]),
                c("Somatic (Pass)", "Germline: both",
                  "Germline: Binom only", "Germline: Betabinom only"))

d$Verdict <- verdict_map[d$Binom_Verdict]
d <- d[!is.na(d$Verdict), ]
d$Verdict <- factor(d$Verdict, levels = names(pal))

theme_clean <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        plot.title       = element_text(face = "bold", size = 13))

## Figure 1 – Distribution of log10(q-value) by verdict group
sub1 <- d[is.finite(d$Binom_Germline_qval_log10), ]
p1 <- ggplot(sub1, aes(x = Binom_Germline_qval_log10, fill = Verdict, colour = Verdict)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(limits = c(-10, 0), oob = squish) +
  labs(title    = paste0("Binomial Germline q-value Distribution  |  Group: ", group_id),
       x        = expression(log[10](q-value)),
       y        = "Density",
       subtitle = "More negative = stronger somatic signal (VAF well below 0.5); values near 0 are consistent with germline VAF") +
  theme_clean +
  guides(fill   = guide_legend(nrow = 2), colour = guide_legend(nrow = 2))

## Figure 2 – Distribution of Rho (overdispersion) by verdict group
sub2 <- d[!is.na(d$Binom_Rho), ]
p2 <- ggplot(sub2, aes(x = Binom_Rho, fill = Verdict, colour = Verdict)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title    = paste0("Betabinomial Rho (Overdispersion) Distribution  |  Group: ", group_id),
       x        = "Rho", y = "Density",
       subtitle = "Higher Rho = greater read-count overdispersion relative to binomial expectation") +
  theme_clean +
  guides(fill   = guide_legend(nrow = 2), colour = guide_legend(nrow = 2))

## Figure 3 – 2D bin density: log10(q-value) vs Rho
sub3 <- d[is.finite(d$Binom_Germline_qval_log10) & !is.na(d$Binom_Rho), ]
sub3 <- sub3[sub3$Binom_Germline_qval_log10 >= -10, ]
p3 <- ggplot(sub3, aes(x = Binom_Germline_qval_log10, y = Binom_Rho)) +
  geom_bin2d(bins = 80) +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",
                         name = "Count", trans = "log10",
                         labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  scale_x_continuous(limits = c(-10, 0), breaks = seq(-10, 0, 2)) +
  scale_y_continuous(limits = c(0, 1),  breaks = seq(0, 1, 0.1)) +
  geom_vline(xintercept = 0,   linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  annotate("text", x = -0.3, y = 0.12, label = "Rho = 0.1", colour = "grey50", size = 3, hjust = 1) +
  labs(title    = paste0("2D Density: q-value vs Rho  |  Group: ", group_id),
       x        = expression(log[10](q-value) ~ " [shown -10 to 0]"), y = "Rho",
       subtitle = "Variants with q-value ~0 and low Rho cluster as somatic candidates") +
  theme_clean +
  theme(legend.position = "right")

pdf(output_pdf, width = 12, height = 7)
print(p1)
print(p2)
print(p3)
dev.off()
cat("Written:", output_pdf, "\n")
