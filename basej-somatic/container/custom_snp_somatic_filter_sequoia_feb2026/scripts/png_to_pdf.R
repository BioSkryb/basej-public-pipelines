args    <- commandArgs(trailingOnly = TRUE)
png_in  <- args[1]
pdf_out <- args[2]
img <- png::readPNG(png_in)
pdf(pdf_out, width = 12, height = 8)
par(mar = c(0, 0, 0, 0))
plot.new()
rasterImage(img, 0, 0, 1, 1)
invisible(dev.off())
