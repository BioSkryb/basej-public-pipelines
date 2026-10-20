library(ohchibi)
library(dplyr)
library(ggtree)
library(optparse)

set.seed(130816)

bskb_col <- c("#12284C", "#1082A2", "#A0CC2C", "#DD14D3", "#F45D34", "#777776", "#FFFFFF")

# Command line options
option_list <- list(
  make_option(c("--metrics_file"), type = "character", default = NULL,
              help = "Path to metrics file (TSV)", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$metrics_file)) {
  stop("Please provide --metrics_file argument.")
}

df <- read.table(opt$metrics_file, header = TRUE, sep = "\t")

p0 <- ggplot(data = df, aes(MAPD_CNV_Log2, SKEW_CNV)) +
  annotate("rect", fill = "#F2B84B", alpha = 0.3, 
           xmin = 0.01, xmax = 0.2,
           ymin = 0, ymax = 0.2)  +
  annotate("rect", fill = "#D9644A", alpha = 0.3, 
           xmin = 0.01, xmax = 0.25,
           ymin = 0.2, ymax = 0.25)  +
  annotate("rect", fill = "#D9644A", alpha = 0.3, 
           xmin = 0.2, xmax = 0.25,
           ymin = 0, ymax = 0.2)  +
  annotate("rect", fill = "#59323C", color = NA, alpha = 0.3, 
           xmin = 0.25, xmax = 1,
           ymin = 0, ymax = 0.2)  +
  annotate("rect", fill = "#5E7348", color = NA, alpha = 0.3, 
           xmin = 0.01, xmax = 0.2,
           ymin = 0.25, ymax = 1)  +
  annotate("rect", fill = "#22221F", color = NA, alpha = 0.1, 
           xmin = 0.2, xmax = 1,
           ymin = 0.25, ymax = 1)  +
  annotate("rect", fill = "#22221F", color = NA, alpha = 0.1, 
           xmin = 0.25, xmax = 1,
           ymin = 0.2, ymax = 0.25)  +
  annotate("text", color = "black",
           x = 0.8, y = 0.1, label = "Noisier", size = 3) +
  annotate("text", color = "black",
           x = 0.1, y = 0.8, label = "Uneven", size = 3) +
  annotate("text", color = "black",
           x = 0.8, y = 0.8, label = "Noisier and Uneven", size = 3) +
  annotate("text", color = "black",
           x = 0.085, y = 0.175, label = "Ideal", size = 3) +
  annotate("text", color = "black",
           x = 0.1, y = 0.225, label = "Border", size = 3) +
  geom_point(fill = "#1082A2", size = 2, shape = 21) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey"),
    axis.ticks.y = element_line(size = unit(0.1, "line")),
    axis.ticks.x = element_line(size = unit(0.1, "line")),
    axis.text.x = element_text(size = 9, color = "grey30", angle = 90, vjust = 0.5, hjust = 0.1),
    axis.text.y = element_text(size = 9, color = "grey30"),
    panel.background = element_blank(),
    axis.line = element_line(color = 'black', size = 0.3),
    panel.spacing.y = unit(0.1, "lines"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9, color = "grey30"),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9)
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0.001, 0.01), oob = squish) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0.01, 0.005)) +
  xlab(label = "MAPD of CNV bins (BioSkryb)") +
  ylab(label = "Unevenness of CNV bins across segments (BioSkryb)")

oh.save.pdf(p = p0, outname = "CNV-Quadrants.pdf", outdir = getwd(), width = 8, height = 8)
ggsave(filename = "CNV-Quadrants_mqc.jpg", plot = p0, path = "./", width = 16, height = 8, units = "in", dpi = 300)