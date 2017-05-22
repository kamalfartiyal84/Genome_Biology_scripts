
##### The script below process HiCUP aligned HiC bam files
##### Generate normalized HiC matrix from bam files
##### Generate normalized density plots

### load R libraries
library(GenomicAlignments)
library(diffHic) 
require(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)

### set location to hicup processed HiC bam files
all_bams <- list.files(pattern="*hicup.bam")

## generate HindIII cut genome fragments
## Make sure to use the same genome version as used for alignment during HiCUP step
hs.frag <- cutGenome("./GRCh37_g1kp2/fasta/GRCh37_g1kp2.fa", "AAGCTT", 4)
hs.param <- pairParam(hs.frag)

### generate h5 file for each bam file
diagnostics <- preparePairs(all_bams[i], hs.param,file="SampleName.h5", dedup=TRUE, minq=20)

### function to compute normalized HiC matrix
hic_norm <- function(input,param,bin.size)
{
   ### binsize 500kb
  message("computing interaction counts....")
  data <- squareCounts(input, hs.param, width=bin.size, filter=1)
  
  # Each row of the count matrix represents the counts for an interaction, 
  # while each column represents a library
  
  message("remove unwanted chromosomes....")
  all_chrs <- c(paste0(rep("chr",22),1:22),"chrX","chrY")
  reg <- regions(data)
  cm <- inflate(data, reg, reg)
  cm_mat <- as.matrix(as.matrix(cm))
  cm_mat[is.na(cm_mat)] <- 0
  ## to remove unwanted chromosomes
  torm1 <- which(!seqnames(reg) %in% all_chrs)
  ## to remove regions with no interactions
  torm2 <- which(rowSums(cm_mat)==0)
  torm <- sort(unique(c(torm1,torm2)))
  cm_mat <- cm_mat[-c(torm),]
  cm_mat <- cm_mat[,-c(torm)]
  reg <- reg[-c(torm)]
  
  message("compute hindIII count in each genome window of interaction....")
  hindIII_count <- matrix(0,length(reg),length(reg))
  row_ind <- matrix(reg$nfrags,length(reg),length(reg),byrow=FALSE)
  col_ind <- matrix(reg$nfrags,length(reg),length(reg),byrow=TRUE)
  hindIII_count <- row_ind+col_ind
  colnames(cm_mat) <- as.character(seqnames(reg))
  rownames(cm_mat) <- as.character(seqnames(reg))
  
  message("normalize interaction data with hindIII count...")
  cm_mat_norm <- cm_mat/hindIII_count
  cm_mat_norm <- round(cm_mat_norm,3)
  colnames(cm_mat_norm) <- paste0(seqnames(reg),":",start(reg),"-",end(reg))
  rownames(cm_mat_norm) <- paste0(seqnames(reg),":",start(reg),"-",end(reg))
  y <- rev(rownames(cm_mat_norm))
  cm_mat_norm_rev <- cm_mat_norm[y,]
  return(cm_mat_norm_rev)
}

SampleName_matnorm <- hic_norm(input="SampleName.h5",param=hs.param,bin.size=5e5)


##### The commands below plots normalized density
library(gplots)
library(marray)
library(RColorBrewer)
cols <- maPalette(low="white", mid="beige", high="blue", k=128)

combMat <- as.matrix(SampleName_matnorm)
densnorm <- quantile(combMat,.98,na.rm=TRUE)
tnorm <- densnorm
combMat[combMat > tnorm] <- tnorm
combMat <- combMat / tnorm
combMat <- round(combMat,3)
reg_chrs <- table(as.character(seqnames(reg)))
ind <- match(all_chrs,names(reg_chrs))
reg_chrs <- reg_chrs[ind]
colSep <- c(0, cumsum(reg_chrs))
rowSep <- c(0, cumsum(rev(reg_chrs)))

png(filename="SampleName_densityPlot.png",width = 9, height = 9, units = 'in', res = 300)
heatmap.2(combMat, Rowv=FALSE, Colv=FALSE,scale="none",col=cols,sepwidth=c(0.0001,0.0001),
          density.info="none",trace="none",labRow="",labCol="",colsep=colSep,sepcolor='grey',rowsep=rowSep)
dev.off()





