args <- commandArgs(trailing = TRUE)
if (length(args) < 4) {
  cat("Usage: plot_copy_number.R sample copy_number_file segments_file pdf_file\n")
  stop("Error: insufficient arguments")
}

library(copynumber)
library(dplyr)

sample <- args[1]
copy_number_file <- args[2]
segments_file <- args[3]
pdf_file <- args[4]

sample <- "A019"
copy_number_file <- list.files(pattern="*.copyNumber.txt")
segments_file <- list.files(pattern="*segmentedCopyNumber.txt")
pdf_file <- "A019.binSize100.pdf"
i=6
data <- read.delim(copy_number_file[i], stringsAsFactors = FALSE)

data <- data %>%
  mutate(pos = (start + end - 1) / 2) %>%
  select(chrom, pos, copynumber)
colnames(data) <- c("chrom", "pos", sample)

segments <- read.delim(segments_file[i], stringsAsFactors = FALSE)
colnames(segments)[5] <- sample
segments$sampleID <- sample

pdf(pdf_file)

plotGenome(data = data, segments = segments,
           sample = 1,
           connect = FALSE,
           seg.lwd = 1.5, cex = 1,
           ylim = c(-6, 6),
           xlab = "chromosomes",
           ylab = expression("log"[2]*" ratio"))

plotSample(data = data, segments = segments,
           layout = c(4, 2),
           sample = 1,
           connect = TRUE,
           seg.lwd = 1.5, cex = 1,
           ylab = expression("log"[2]*" ratio"))

dev.off()
