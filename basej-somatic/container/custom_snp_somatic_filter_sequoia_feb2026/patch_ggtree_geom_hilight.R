# ggtree 3.6.2 + ggplot2 3.5.x: ggplot2 removed check_linewidth from its namespace.
# Inline the helper used in current YuLab-SMU/ggtree (see geom_hilight.R in ggtree >= 3.12).
f <- "/tmp/ggtree/R/geom_hilight.R"
lines <- readLines(f)
pat <- "check_linewidth <- getFromNamespace('check_linewidth', 'ggplot2')"
i <- which(lines == pat)
if (length(i) != 1L) stop("patch_ggtree: expected exactly one line matching: ", pat)
insert <- c(
  ".check_linewidth <- function(data, name) {",
  "  if (is.null(data$linewidth) && !is.null(data$size)) {",
  "    warning(paste0(",
  "      \"Using the `size` aesthetic with \", name, \" was deprecated in ggplot2 3.4.0.\\n\",",
  "      \"Please use the `linewidth` aesthetic instead.\"",
  "    ))",
  "    data$linewidth <- data$size",
  "  }",
  "  data",
  "}"
)
lines <- c(lines[seq_len(i - 1L)], insert, lines[seq(i + 1L, length(lines))])
lines <- gsub("check_linewidth(data,", ".check_linewidth(data,", lines, fixed = TRUE)
writeLines(lines, f)
