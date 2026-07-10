# theme_ohchibi_pubr.R — ggplot2 theme tuned for publication-style figures
# (dotted major grid, grey30 tick labels, bottom legend).
# Axis titles, legend title/text, strip (facet) text, plot title/subtitle, and axis tick
# labels: size = 9 (absolute points).
#
# Usage:
#   library(ggplot2)
#   source("/home/ubuntu/projects/theme_ohchibi_pubr.R")
#   ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_ohchibi_pubr()

theme_ohchibi_pubr <- function(base_size = 13, base_family = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for theme_ohchibi_pubr()", call. = FALSE)
  }
  half_line <- base_size / 2
  base <- ggplot2::theme_bw(base_size = base_size, base_family = base_family)
  overlay <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.grid.major.x = ggplot2::element_line(linetype = "dotted", colour = "grey"),
    panel.grid.major.y = ggplot2::element_line(linetype = "dotted", colour = "grey"),
    panel.grid.minor.y = ggplot2::element_line(colour = "grey93", linewidth = 0.2),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.3),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
    axis.text.x = ggplot2::element_text(
      colour = "grey30", face = "plain", size = 9,
      angle = 40, hjust = 1, vjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      colour = "grey30", face = "plain", size = 9
    ),
    axis.title.x = ggplot2::element_text(
      colour = "black", face = "plain", size = 9
    ),
    axis.title.y = ggplot2::element_text(
      colour = "black", face = "plain", size = 9,
      angle = 90, vjust = 1
    ),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "plain", size = 9, margin = ggplot2::margin(b = 4)
    ),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_text(
      face = "plain", colour = "black", size = 9
    ),
    legend.text = ggplot2::element_text(
      colour = "black", face = "plain", size = 9
    ),
    legend.key.size = grid::unit(0.70, "cm"),
    legend.key.width = grid::unit(0.70, "cm"),
    legend.spacing.x = grid::unit(0.4, "cm"),
    legend.margin = ggplot2::margin(t = 8, b = 4),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "plain", size = 9),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "plain", size = 9),
    plot.margin = ggplot2::margin(half_line, half_line, half_line, half_line)
  )
  ggplot2::`%+replace%`(base, overlay)
}
