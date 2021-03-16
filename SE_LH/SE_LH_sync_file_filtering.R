## 
## This script analyzes output from Popoolation2. 
##

library(poolfstat)
library(rtracklayer)
library(scales)
library(Rfast)

## set population colors
se.col <- '#992913'
lh.col <- '#D6A929'

setwd('/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/snps_and_gene_expression/variant_calling/pop_pools/20_05_02_SE_LH_10_pools/')

###### 1. READ IN DATA #####--------------------------------------------------------------------
## read in popoolation2 sync file
## man page says poolsizes should be the haploid pool sizes, so 40?
sync.data <- popsync2pooldata(sync.file='/Users/Avril/Desktop/.scratch/pop_comp/popoolation2_10_pools/no_zeroes_all_chrom.sync', poolsizes=c(8,8,8,8,8,8,8,8,8,8), noindel=TRUE)
## using 'no_zeroes' file = '0:0:0:0:0:0' not allowed for ANY sample (i.e., at least 1X covg of at least 1 allele in all samples)
## processing no_zeroes_all_chrom.sync = ~51 minutes,  271,692 SNPs/48 million lines

# read in annotation information for later annotation of SNPs
full.assembly <- readGFF("/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/stringtie_featureCounts_deseq2_DGE_analyses/stringtie_all_merged.annotated.gtf")
full.assembly$full.name <- paste0(full.assembly$gene_name,'|',full.assembly$transcript_id,'|',full.assembly$xloc)

# read in gene list from published assembly
descrip <- read.csv('/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/stringtie_featureCounts_deseq2_DGE_analyses/ICSASG_v2_gene_list.csv')

sync.data@poolnames <- c('SE.16','SE.17','SE.18','SE.19','SE.20','LH.13','LH.14','LH.15','LH.21','LH.22')

##### 2. Set required coverage for SNPs #####--------------------------------------------------------------------
covg <- sync.data@readcoverage
rownames(covg) <- paste0(sync.data@snp.info[,1],'|',sync.data@snp.info[,2])
colnames(covg) <- c('SE.16','SE.17','SE.18','SE.19','SE.20','LH.13','LH.14','LH.15','LH.21','LH.22')

##### !!! #####
## set required coverage
## require that 4/5 pools per population have nX covg.
n <- 8 ## at 8X for 4/5 samples --> 162,204 (SE) & 160,397 (LH)

se.covg <- covg[,c(1:5)]
lh.covg <- covg[,c(6:10)]

dim(se.covg[rowSums(se.covg >= n) >= 4,])
dim(lh.covg[rowSums(lh.covg >= n) >= 4,])

se.covg <- se.covg[rowSums(se.covg >= n) >= 4,]
lh.covg <- lh.covg[rowSums(lh.covg >= n) >= 4,]

## ID SNPs covered in both SE and LH
keep.snps <- rownames(se.covg[which(rownames(se.covg) %in% rownames(lh.covg)),])
req.covg <- covg[which(rownames(covg) %in% keep.snps),] ## 135,839 SNPs at 8X for 4/5 samples * 2 pops

## check variance in read depth across samples for all loci
max(rowVars(req.covg, std=TRUE))
min(rowVars(req.covg, std=TRUE))
hist(rowVars(req.covg, std=TRUE), xlim=c(0,50), breaks=100)
hist(rowVars(req.covg, std=TRUE)/rowMeans(req.covg)) ## plot CoV hist
req.covg[which(rowVars(req.covg, std=TRUE)/rowMeans(req.covg) > 3.0),] ## no SNPs with CoV > 3 SD from the mean

##### QUANTIFYING DEPTH VARIATION IN ALL SNPS, FILTERED SNPS, AND OUTLIER SNPS #####
head(req.covg)
## compile 10-pool data to be filtered
pool.10.data <- as.data.frame(sync.data@snp.info)
pool.10.data <- cbind(pool.10.data, sync.data@refallele.readcount)
pool.10.data <- cbind(pool.10.data, sync.data@readcoverage)
colnames(pool.10.data) <- c('chrom','pos','pool.ref','pool.alt','ref.SE.16','ref.SE.17','ref.SE.18','ref.SE.19','ref.SE.20','ref.LH.13','ref.LH.14','ref.LH.15','ref.LH.21','ref.LH.22','dp.SE.16','dp.SE.17','dp.SE.18','dp.SE.19','dp.SE.20','dp.LH.13','dp.LH.14','dp.LH.15','dp.LH.21','dp.LH.22')
pool.10.data$loc <- paste0(pool.10.data$chrom,'|',pool.10.data$pos)
pool.10.data <- pool.10.data[which(pool.10.data$loc %in% rownames(req.covg)),]
pool.10.data <- pool.10.data[pool.10.data$chrom %in% unique(pool.10.data$chrom)[2:30],] ## keep only SNPs located on assembled nuclear chromosomes
write.table(pool.10.data, '/Users/Avril/Desktop/10_pool_109656_SNPs_allele_counts.txt',
  sep='\t', row.names=FALSE, quote=FALSE) ## this file is read into 3 population_analyses.R

filt.covg <- req.covg
## Calculate CoV for read depth using all 10 pools combined (results in a more conservative cut-off than calculating for each population separately)
filt.stats <- cbind(rowMeans(filt.covg), rowVars(filt.covg, std=TRUE)/rowMeans(filt.covg))
colnames(filt.stats) <- c('rowmean','rowCoV')
rownames(filt.stats) <- rownames(filt.covg)
stdev <- sd(filt.stats[,2])
lowlim <- mean(filt.stats[,2])-(3*stdev)
uplim <- mean(filt.stats[,2])+(3*stdev) ## upper limit for removing CoV outliers = mean CoV + 3 SD

## plot CoV for all SNPs
all <- 'grey40'
out <- '#FFBF00'
# pdf('/Users/Avril/Desktop/CoV_filtering.pdf', width=6, height=5)
plot(filt.stats[,1], filt.stats[,2], pch=19, col=alpha(all), cex=0.75,
     xlab='Mean SNP read depth', ylab='Coefficient of variation', main='10 pools')
  points(knn.stats[,1], knn.stats[,2], pch=23, col='black', bg=out, cex=0.75)
  legend('topright', inset=0.05, legend=c('Filtered SNPs','Outlier SNPs'), col=c(all,'black'), pt.bg=c(NA,out), pch=c(19, 23), pt.cex=0.75)
  lines(x=c(min(filt.stats[,1]), max(filt.stats[,1])), y=c(uplim, uplim), col='black', lty=2)
#dev.off()
