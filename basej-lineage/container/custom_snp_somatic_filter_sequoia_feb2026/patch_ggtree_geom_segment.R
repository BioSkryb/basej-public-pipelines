# ggtree 3.6.2 + ggplot2 3.5.x: ggplot2's empty() now relies on is_waiver(); getFromNamespace("empty","ggplot2")
# breaks when that runs in ggtree's environment. Match YuLab-SMU/ggtree master R/geom_segment.R.
f <- "/tmp/ggtree/R/geom_segment.R"
text <- paste(readLines(f), collapse = "\n")
old <- paste0(
  "empty <- getFromNamespace(\"empty\", \"ggplot2\")\n",
  "`%||%` <- getFromNamespace(\"%||%\", \"ggplot2\")"
)
if (!grepl(old, text, fixed = TRUE)) stop("patch_geom_segment: expected tail not found")
new <- paste0(
  "is_waiver <- function(x) {\n",
  "  inherits(x, \"waiver\")\n",
  "}\n",
  "\n",
  "empty <- function(df) {\n",
  "  is.null(df) || nrow(df) == 0 || ncol(df) == 0 || is_waiver(df)\n",
  "}\n",
  "\n",
  "`%||%` <- function(a, b) {\n",
  "  if (!is.null(a)) a else b\n",
  "}"
)
text <- sub(old, new, text, fixed = TRUE)
writeLines(strsplit(text, "\n", fixed = TRUE)[[1]], f)
