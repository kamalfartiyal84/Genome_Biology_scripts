args <- commandArgs(trailing = TRUE)
if (length(args) < 2) {
  cat("Usage: apply_restriction_site_correction.R Rdata_file restriction_site_file corrected_copy_number_file\n")
  stop("Error: insufficient arguments")
}
rdata_file <- args[1]
restriction_site_file <- args[2]
corrected_copy_number_file <- args[3]

library(Biobase)
library(QDNAseq)
library(dplyr)

load(rdata_file)

counts <- assayDataElement(readCounts, "counts")

restriction_sites <- read.delim(restriction_site_file, stringsAsFactors = FALSE)
mean_restriction_site_count <- restriction_sites %>% filter(count >= 10 && count <= 50) %>% select(count) %>% unlist %>% mean
restriction_sites <- restriction_sites %>% mutate(scale_factor = mean_restriction_site_count / ifelse(count == 0, mean_restriction_site_count, count))
excluded <- restriction_sites %>% filter(count < 10 || count > 50)
fData(readCounts)[excluded$location, "use"] <- FALSE

restriction_sites <- left_join(data_frame(location = rownames(counts)), restriction_sites, by = "location")

assayDataElement(readCounts, "counts") <- counts * restriction_sites$scale_factor

readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
readCountsCorrected <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsCorrected)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

copyNumbers <- data.frame(
  location = rownames(copyNumbersSegmented),
  chrom = featureData(copyNumbersSegmented)$chromosome,
  start = as.integer(featureData(copyNumbersSegmented)$start),
  end = as.integer(featureData(copyNumbersSegmented)$end),
  copynumber = round(log2(assayData(copyNumbersSegmented)$copynumber)[,1], digits = 3),
  segmented = round(log2(assayData(copyNumbersSegmented)$segmented)[,1], digits = 3)
)

options(scipen = 999)
write.table(copyNumbers, file = corrected_copy_number_file, row.names = FALSE, quote = FALSE, sep = "\t")
