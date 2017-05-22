library(Biobase)
library(QDNAseq)
library(CGHcall)

args <- commandArgs(trailing = TRUE)
bamfile <- args[1]
id <- args[2]
name <- args[3]
prefix <- args[4]
binSize <- as.numeric(args[5])

bins <- readRDS(paste("./QDNAseq/binAnnotations", 100, ".rds", sep = ""))

readCounts <- binReadCounts(bins, bamfiles = bamfile)
exportBins(readCounts, file = paste(prefix, "readCounts.txt", sep = "."), logTransform = FALSE)

png(paste(prefix, "readCounts.png", sep = "."), width = 1200, height = 800)
plot(readCounts, logTransform = FALSE, ylim = c(-10, binSize * 15), main = paste(id, name))
highlightFilters(readCounts, logTransform = FALSE, residual = TRUE, blacklist = TRUE)
dev.off()

readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)


png(paste(prefix, "isobar.png", sep = "."), width = 1200, height = 800)
isobarPlot(readCountsFiltered)
dev.off()

readCountsCorrected <- estimateCorrection(readCountsFiltered)

rawReadCounts <- data.frame(
	location = rownames(readCountsCorrected),
	chrom = featureData(readCountsCorrected)$chromosome,
	start = as.integer(featureData(readCountsCorrected)$start),
	end = as.integer(featureData(readCountsCorrected)$end),
	count = assayData(readCountsCorrected)$counts[,1]
)
write.table(rawReadCounts, paste(prefix, "rawReadCounts.txt", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

fittedReadCounts <- data.frame(
	location = rownames(readCountsCorrected),
	chrom = featureData(readCountsCorrected)$chromosome,
	start = as.integer(featureData(readCountsCorrected)$start),
	end = as.integer(featureData(readCountsCorrected)$end),
	count = round(assayData(readCountsCorrected)$fit[,1], digits = 3)
)
write.table(fittedReadCounts, paste(prefix, "fittedReadCounts.txt", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

png(paste(prefix, "noise.png", sep = "."), width = 1200, height = 800)
noisePlot(readCountsCorrected)
dev.off()

copyNumbers <- correctBins(readCountsCorrected)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

exportBins(copyNumbersSmooth, file = paste(prefix, "igv", sep = "."), format = "igv")

copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

png(paste(prefix, "copyNumberSegmented.png", sep = "."), width = 1200, height = 800)
plot(copyNumbersSegmented, main = paste(id, name))
dev.off()


copyNumbers <- data.frame(
	location = rownames(copyNumbersSegmented),
	chrom = featureData(copyNumbersSegmented)$chromosome,
	start = as.integer(featureData(copyNumbersSegmented)$start),
	end = as.integer(featureData(copyNumbersSegmented)$end),
	copynumber = round(log2(assayData(copyNumbersSegmented)$copynumber)[,1], digits = 3),
	segmented = round(log2(assayData(copyNumbersSegmented)$segmented)[,1], digits = 3)
)

options(scipen = 999)
write.table(copyNumbers, file = paste(prefix, "copyNumber.txt", sep = "."), row.names = FALSE, quote = FALSE, sep = "\t")

copyNumbersCalled <- callBins(copyNumbersSegmented, nclass = 3)
cgh <- makeCgh(copyNumbersCalled)

png(paste(prefix, "copyNumberCalls.png", sep = "."), width = 1200, height = 800)
plot(copyNumbersCalled, main = paste(id, name))
dev.off()

copyNumberCalls <- data.frame(
	location = rownames(cgh),
	chrom = featureData(cgh)$Chromosome,
	start = as.integer(featureData(cgh)$Start),
	end = as.integer(featureData(cgh)$End),
	copynumber = round(copynumber(cgh)[,1], digits = 3),
	segmented = round(segmented(cgh)[,1], digits = 3),
	call = calls(cgh)[,1],
	probnorm = round(probnorm(cgh)[,1], digits = 2),
	probgain = round(probgain(cgh)[,1], digits = 2),
	probloss = round(probloss(cgh)[,1], digits = 2)
)

options(scipen = 999)
write.table(copyNumberCalls, file = paste(prefix, "copyNumberCalls.txt", sep = "."), row.names = FALSE, quote = FALSE, sep = "\t")

save.image(file = paste(prefix, "Rdata", sep = "."))

sessionInfo()

