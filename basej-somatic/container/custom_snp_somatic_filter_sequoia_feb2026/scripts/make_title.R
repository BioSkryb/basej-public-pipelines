args    <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
title   <- args[2]
group   <- args[3]
pdf(outfile, width = 11, height = 8.5)
par(mar = c(0, 0, 0, 0))
plot.new()
rect(0, 0, 1, 1, col = "#003366", border = NA)
text(0.5, 0.55, title, cex = 2.6, font = 2, col = "white",   adj = c(0.5, 0.5))
text(0.5, 0.35, group, cex = 1.4, font = 3, col = "#aaccff", adj = c(0.5, 0.5))
invisible(dev.off())
