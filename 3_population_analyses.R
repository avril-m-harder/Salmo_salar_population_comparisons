setwd('/Users/Avril/Desktop/all_ch_3_dirs/02_output_tables_all_analyses/')

library(rtracklayer)
library(vcfR)
library(scales)
library(ggrepel)
library(dplyr)
library(ade4)
library(adegenet)
library(ape)
library(dendextend)
library(polysat)

`%notin%` <- Negate(`%in%`)

## set comparison colors
lc.vs.se <- '#47304a'
lc.vs.lh <- '#229284'  
se.vs.lh <- '#e27b0f'
## set population colors
lc.col <- '#1571a9'
se.col <- '#992913'
lh.col <- '#D6A929'

##### process Lake Champlain data #####
## read in VCF and gtypes generated in LC_vcf_filtering.R based on 
## file created by variant_filtration.sh script
vcf.filt.snps <- read.vcfR('/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/snps_and_gene_expression/variant_calling/20_12_15_testing_scripts_for_ms_upload/post_filtering_36_samples.vcf')
lc.gtypes <- read.csv('/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/snps_and_gene_expression/variant_calling/20_12_15_testing_scripts_for_ms_upload/20_12_28_genotype_data_36_samples.csv', row.names=1, check.names=FALSE)

## calculate LC genotype frequencies to get allele frequencies
OUT <- NULL
for(i in 1:nrow(lc.gtypes)){
  if(length(which(is.na(lc.gtypes[i,]))) < 16){ ## if <= 16 samples are missing genotypes (require 20 genotypes),
    ## calculate number of each genotype
    num.hom.ref <- length(lc.gtypes[i,which(lc.gtypes[i,]=='0')])
    num.het <- length(lc.gtypes[i,which(lc.gtypes[i,]=='1')])
    num.hom.alt <- length(lc.gtypes[i,which(lc.gtypes[i,]=='2')])
    tot <- sum(num.hom.ref, num.het, num.hom.alt)
    ## calculate allele frequencies from genotypes
    ref.af <- (((num.hom.ref*2) + (num.het)) / (tot*2)) ## number of ref allele in hom.ref + number in het / total # of alleles in indivs w/ genotypes
    alt.af <- (((num.hom.alt*2) + (num.het)) / (tot*2)) ## same for alt allele
    ## and store
    save <- c(rownames(lc.gtypes)[i], tot, ref.af, alt.af)
    OUT <- rbind(OUT, save)
  }
  print(i/nrow(lc.gtypes))
}
dim(OUT) 

lc.afs <- as.data.frame(OUT)
rownames(lc.afs) <- lc.afs$V1
lc.afs$chrom <- do.call(rbind, strsplit(rownames(lc.afs), split='|', fixed=TRUE))[,1]
lc.afs$pos <- do.call(rbind, strsplit(rownames(lc.afs), split='|', fixed=TRUE))[,2]
lc.afs <- lc.afs[,-1]
colnames(lc.afs) <- c('samps','lc.ref.af','lc.alt.af','chrom','pos')

## get allele depth and overall depth information
gt.info <- vcf.filt.snps@gt
loci <- vcf.filt.snps@fix[,c(1,2)]
loc <- paste0(loci[,1],'|',loci[,2])
rownames(gt.info) <- loc
gt.info <- gt.info[which(rownames(gt.info) %in% rownames(lc.afs)),]

##### Prepping for SE/LH comparisons #####
## read in SE/LH data just to limit down LC loci to shared loci
pools.10.dp <- read.table('10_pool_109656_SNPs_allele_counts.txt', sep='\t', header=TRUE)
gt.info <- gt.info[which(rownames(gt.info) %in% pools.10.dp$loc),]
## get rid of FORMAT column
gt.info <- gt.info[,-1]
## reformat colnames
samps <- colnames(gt.info)
split.samps <- do.call(rbind, strsplit(samps, split='_', fixed=TRUE))[,4]
colnames(gt.info) <- do.call(rbind, strsplit(split.samps, split='.', fixed=TRUE))[,1]
## loop over columns in gt.info, saving allele and total depth information
REF.DP <- NULL
TOTAL.DP <- NULL
for(i in 1:ncol(gt.info)){
  ads <- do.call(rbind, strsplit(gt.info[,i], split=':', fixed=TRUE))[,2]
  tot <- as.numeric(do.call(rbind, strsplit(gt.info[,i], split=':', fixed=TRUE))[,3])
  ref.dp <- as.numeric(do.call(rbind, strsplit(ads, split=',', fixed=TRUE))[,1])
  REF.DP <- cbind(REF.DP, ref.dp)
  TOTAL.DP <- cbind(TOTAL.DP, tot)
}
rm(ref.dp)
colnames(REF.DP) <- colnames(gt.info)
colnames(TOTAL.DP) <- colnames(gt.info)
ref.dp <- REF.DP
rm(REF.DP)
total.dp <- TOTAL.DP
rm(TOTAL.DP)
rownames(ref.dp) <- rownames(gt.info) 
rownames(total.dp) <- row.names(gt.info)
## get reference allele identity information (for matching with SE/LH alleles)
allele.ids <- as.data.frame(vcf.filt.snps@fix[,c(1,2,4,5)])
allele.ids$loc <- paste0(allele.ids$CHROM,'|',allele.ids$POS)
allele.ids <- allele.ids[which(allele.ids$loc %in% rownames(ref.dp)),]


##### process SE/LH data #####
## get REF and total depth information for 10 pools - written from SE_LH_only_10_pools.R for all
## 109,656 SNPs meeting pool depth filtering criteria
pools.10.dp <- read.table('10_pool_109656_SNPs_allele_counts.txt', sep='\t', header=TRUE) ## allele identities need to be fixed

##### calculate 10-pool REF AFs and correct based on ref allele identity #####
## Correcting reference allele ID is necessary because poolfstat, used to filter SE/LH data,
## artbitrarily selects an allele as the "reference" allele
pools.10.dp <- pools.10.dp[which(pools.10.dp$loc %in% allele.ids$loc),]
allele.ids <- allele.ids[which(allele.ids$loc %in% pools.10.dp$loc),]
allele.ids <- allele.ids[order(allele.ids$loc),]
pools.10.dp <- pools.10.dp[order(pools.10.dp$loc),]
all(allele.ids$loc == pools.10.dp$loc) ## check that locus order matches
pools.10.dp <- cbind(pools.10.dp, allele.ids[,c(3,4)])
colnames(pools.10.dp)[c(26,27)] <- c('lc.ref','lc.alt')

## correct ref allele identity for SE/LH pools
MATCH <- as.character(pools.10.dp$lc.ref) == as.character(pools.10.dp$pool.ref)
table(MATCH) ## allele identities do not all match between LC and SE/LH
pools.10.dp <- cbind(pools.10.dp, MATCH)
## pull out loci that need to be corrected
fix.allele.IDs <- pools.10.dp[which(pools.10.dp$MATCH == FALSE),]
fix.allele.IDs$loc <- as.character(fix.allele.IDs$loc)
## check to make sure that ALT allele in SE/LH matches the LC REF allele
ALT.MATCH <- as.character(fix.allele.IDs$lc.ref) == as.character(fix.allele.IDs$pool.alt)
fix.allele.IDs <- cbind(fix.allele.IDs, ALT.MATCH)
table(fix.allele.IDs$ALT.MATCH) ## 11 SNPs have 2 distinct sets of alleles between the two data sets (should discard dowstream,
## flag here)
## correct snp.fsts REF AF information in fix.allele.IDs, then replace information in pools.10.dp
for(i in 1:nrow(fix.allele.IDs)){ ## for every SNP where REF and ALT need to be switched for SE/LH
  fix.allele.IDs$ref.SE.16[i] <- (fix.allele.IDs$dp.SE.16[i] - fix.allele.IDs$ref.SE.16[i])
  fix.allele.IDs$ref.SE.17[i] <- (fix.allele.IDs$dp.SE.17[i] - fix.allele.IDs$ref.SE.17[i])
  fix.allele.IDs$ref.SE.18[i] <- (fix.allele.IDs$dp.SE.18[i] - fix.allele.IDs$ref.SE.18[i])
  fix.allele.IDs$ref.SE.19[i] <- (fix.allele.IDs$dp.SE.19[i] - fix.allele.IDs$ref.SE.19[i])
  fix.allele.IDs$ref.SE.20[i] <- (fix.allele.IDs$dp.SE.20[i] - fix.allele.IDs$ref.SE.20[i])
  fix.allele.IDs$ref.LH.13[i] <- (fix.allele.IDs$dp.LH.13[i] - fix.allele.IDs$ref.LH.13[i])
  fix.allele.IDs$ref.LH.14[i] <- (fix.allele.IDs$dp.LH.14[i] - fix.allele.IDs$ref.LH.14[i])
  fix.allele.IDs$ref.LH.15[i] <- (fix.allele.IDs$dp.LH.15[i] - fix.allele.IDs$ref.LH.15[i])
  fix.allele.IDs$ref.LH.21[i] <- (fix.allele.IDs$dp.LH.21[i] - fix.allele.IDs$ref.LH.21[i])
  fix.allele.IDs$ref.LH.22[i] <- (fix.allele.IDs$dp.LH.22[i] - fix.allele.IDs$ref.LH.22[i])
  fix.allele.IDs$AF.FIXED[i] <- TRUE ## and mark it as fixed
  print(i/nrow(fix.allele.IDs))
}
## add back into pools.10.dp
pools.10.dp <- pools.10.dp[which(pools.10.dp$loc %notin% fix.allele.IDs$loc),]
pools.10.dp$AF.FIXED <- NA
odd.snps <- fix.allele.IDs[which(fix.allele.IDs$ALT.MATCH == FALSE), 'loc']
fix.allele.IDs <- fix.allele.IDs[,-29] ## remove ALT.MATCH column (covered by odd.snps, below)
pools.10.dp <- rbind(pools.10.dp, fix.allele.IDs)
## add note for SNPs with 2 distinct allele sets in the two data sets
pools.10.dp[which(pools.10.dp$loc %in% odd.snps), 'FLAGGED'] <- TRUE
## remove flagged SNPs (tri-allelic across all 3 populations)
pools.10.dp <- pools.10.dp[is.na(pools.10.dp$FLAGGED),]

## calculate identity-corrected reference allele frequencies
pools.10.dp$af.SE.16 <- pools.10.dp$ref.SE.16/pools.10.dp$dp.SE.16
pools.10.dp$af.SE.17 <- pools.10.dp$ref.SE.17/pools.10.dp$dp.SE.17
pools.10.dp$af.SE.18 <- pools.10.dp$ref.SE.18/pools.10.dp$dp.SE.18
pools.10.dp$af.SE.19 <- pools.10.dp$ref.SE.19/pools.10.dp$dp.SE.19
pools.10.dp$af.SE.20 <- pools.10.dp$ref.SE.20/pools.10.dp$dp.SE.20
pools.10.dp$af.LH.13 <- pools.10.dp$ref.LH.13/pools.10.dp$dp.LH.13
pools.10.dp$af.LH.14 <- pools.10.dp$ref.LH.14/pools.10.dp$dp.LH.14
pools.10.dp$af.LH.15 <- pools.10.dp$ref.LH.15/pools.10.dp$dp.LH.15
pools.10.dp$af.LH.21 <- pools.10.dp$ref.LH.21/pools.10.dp$dp.LH.21
pools.10.dp$af.LH.22 <- pools.10.dp$ref.LH.22/pools.10.dp$dp.LH.22

## ID overlap in SE/LH/LC filtered loci
locs <- pools.10.dp[which(pools.10.dp$loc %in% rownames(lc.afs)), 'loc'] ## 23,673 SNPs
## limit down to this list of loci
lc.afs <- lc.afs[which(rownames(lc.afs) %in% locs),]
lc.gtypes <- lc.gtypes[which(rownames(lc.gtypes) %in% locs),]
ref.dp <- ref.dp[which(rownames(ref.dp) %in% locs),]
total.dp <- total.dp[which(rownames(total.dp) %in% locs),]
pools.10.dp <- pools.10.dp[which(pools.10.dp$loc %in% locs),]
rm(allele.ids)
rm(fix.allele.IDs)

##### Calculating AF values #####
## combine LC AFs calculated from individual genotypes and SE/LH AFs calculated from combining pool AFs
se.afs <- pools.10.dp[,c(31:35)]
lh.afs <- pools.10.dp[,c(36:40)]
rownames(se.afs) <- pools.10.dp$loc
rownames(lh.afs) <- pools.10.dp$loc
## calculate population AFs from pool AFs for SE/LH
for(i in 1:nrow(se.afs)){
  se.afs$se.ref.af[i] <- ((sum(se.afs[i,c(1:5)])) / (length(which(!is.na(se.afs[i,c(1:5)])))))
}
for(i in 1:nrow(lh.afs)){
  lh.afs$lh.ref.af[i] <- ((sum(lh.afs[i,c(1:5)])) / (length(which(!is.na(lh.afs[i,c(1:5)])))))
}
se.afs <- se.afs[order(rownames(se.afs)),]
lh.afs <- lh.afs[order(rownames(lh.afs)),]
## get LC REF AF calculated from individual genotypes
lc.afs <- lc.afs[order(rownames(lc.afs)),]
all(rownames(lc.afs) == rownames(se.afs))
all(rownames(lc.afs) == rownames(lh.afs))
## combine 3-pop REF AFs
pop.afs <- cbind(as.numeric(paste(lc.afs$lc.ref.af)), se.afs$se.ref.af, lh.afs$lh.ref.af)
pop.afs <- as.data.frame(pop.afs)
rownames(pop.afs) <- rownames(lc.afs)
colnames(pop.afs) <- c('lc.ref.af','se.ref.af','lh.ref.af')
pop.afs$chrom <- do.call(rbind, strsplit(rownames(pop.afs), split='|', fixed=TRUE))[,1]
pop.afs$pos <- do.call(rbind, strsplit(rownames(pop.afs), split='|', fixed=TRUE))[,2]
# write.table(pop.afs, '/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/snps_and_gene_expression/variant_calling/20_12_15_testing_scripts_for_ms_upload/pop_afs_23673_SNPs.txt', sep='\t', quote=FALSE)

##### MAF filtering #####
## calculate MAF
ref.afs <- (pop.afs$lc.ref.af*36 + pop.afs$se.ref.af*20 + pop.afs$lh.ref.af*20) / 76
MAFS <- NULL
for(i in 1:length(ref.afs)){
  if(ref.afs[i] <= 0.5){
    MAFS <- c(MAFS, ref.afs[i])
  }
  else{
    MAFS <- c(MAFS, (1-ref.afs[i]))
  }
}
pop.afs$maf <- MAFS
## only keep loci with MAF > 1/72 (excluding ~singletons)
dim(pop.afs[which(pop.afs$maf > (1/72)),]) ## 23,291 SNPs
pop.afs <- pop.afs[which(pop.afs$maf > (1/72)),]
## limit down to this list of loci
locs <- rownames(pop.afs)
lc.afs <- lc.afs[which(rownames(lc.afs) %in% locs),]
lc.gtypes <- lc.gtypes[which(rownames(lc.gtypes) %in% locs),]
ref.dp <- ref.dp[which(rownames(ref.dp) %in% locs),]
total.dp <- total.dp[which(rownames(total.dp) %in% locs),]
pools.10.dp <- pools.10.dp[which(pools.10.dp$loc %in% locs),]
rm(se.afs)
rm(lh.afs)

##### Generate LC pools to calculate interpopulation genetic distance #####
## generate pools randomly
pools <- sample(1:36, size=36, replace=FALSE)
p.size <- 4
step <- p.size - 1 
## calculate per-pool reference allele frequency and save in one dataframe
lc.pool.afs <- NULL
i <- 1
while(i < 36){
  lc.pool.afs <- cbind(lc.pool.afs, (rowSums(ref.dp[,c(pools[i:(i+step)])]/rowSums(total.dp[,c(pools[i:(i+step)])]))))
  i <- i+p.size
  print(i)
}
colnames(lc.pool.afs) <- c('af.LC.1','af.LC.2','af.LC.3','af.LC.4','af.LC.5','af.LC.6','af.LC.7','af.LC.8','af.LC.9')

lc.pool.afs <- lc.pool.afs[order(rownames(lc.pool.afs)),]
pools.10.dp <- pools.10.dp[order(pools.10.dp$loc),]
all(rownames(lc.pool.afs) == pools.10.dp$loc) ## make sure locus order matches
all.afs <- cbind(pools.10.dp[,c(31:40)], lc.pool.afs)
rownames(all.afs) <- rownames(lc.pool.afs)

##### Format data for 19 pools for Fst calculation with Jacob's script ##### 
## Note: per Jacob's recommendation, all depths will be set equal to the number of individuals
## per population, because for all loci, read depth exceeds the number of individuals sequenced;
## write file with necessary information for fst_from_pooled_freqs.pl
pool.dp <- seq(from=4, to=4, length.out=nrow(all.afs))
fst.input <- cbind(rownames(all.afs), all.afs[,1], pool.dp, all.afs[,2], pool.dp, all.afs[,3], pool.dp, all.afs[,4], pool.dp, all.afs[,5],
                   pool.dp, all.afs[,6], pool.dp, all.afs[,7], pool.dp, all.afs[,8], pool.dp, all.afs[,9], pool.dp, all.afs[,10],
                   pool.dp, all.afs[,11], pool.dp, all.afs[,12], pool.dp, all.afs[,13], pool.dp, all.afs[,14], pool.dp, all.afs[,15],
                   pool.dp, all.afs[,16], pool.dp, all.afs[,17], pool.dp, all.afs[,18], pool.dp, all.afs[,19], pool.dp)
colnames(fst.input) <- c('Allele','Freq1','Size1','Freq2','Size2','Freq3','Size3','Freq4','Size4','Freq5','Size5',
                         'Freq6','Size6','Freq7','Size7','Freq8','Size8','Freq9','Size9','Freq10','Size10',
                         'Freq11','Size11','Freq12','Size12','Freq13','Size13','Freq14','Size14','Freq15','Size15',
                         'Freq16','Size16','Freq17','Size17','Freq18','Size18','Freq19','Size19')
## maybe write in batches of loci to make it work locally?
## write a file of Freq/Size information for each population pair (n=342), then analyze each using a bash loop + Jacob's script
f.cols <- seq(2,38,2) ## vector of column indices for Freq information
for(i in f.cols){
  comps <- f.cols[which(f.cols > i)]
  for(j in comps){
    sub <- fst.input[,c(1,i,i+1,j,j+1)]
    write.table(sub, paste0('/Users/Avril/Desktop/19_pool_analysis/',i/2,'_vs_',j/2,'.txt'), row.names=FALSE, sep='\t', quote=FALSE)
  }
}
# # >>> # run Jacob's script: perl fst_from_pooled_freqs.pl -f allele_freqs.txt -t -1 -d 4
# # <<< # read in output from that script and store average Fst values in matrix
poolnums <- seq(1:19)
OUT <- matrix(nrow=19, ncol=19)
for(i in poolnums){
  comps <- poolnums[which(poolnums > i)]
  print(paste0('i = ',i))
  for(j in comps){
    fsts <- read.table(paste0('/Users/Avril/Desktop/19_pool_analysis/fst_output/High_Fsts_',i,'_vs_',j,'.txt'),
                       header=TRUE, sep='\t', check.names=FALSE)
    fsts[which(fsts[,2] < 0), 2] <- 0 ## convert negative values to 0 or hclust doesn't work
    avg <- mean(fsts[,2])
    OUT[j,i] <- avg ## fill in the matrix
  }
}
## add pool info to row/col names
colnames(OUT) <- colnames(all.afs)
rownames(OUT) <- colnames(all.afs)
## convert to a dist object for plotting
pwise.pool.fst <- as.dist(OUT)
dend <- hclust(pwise.pool.fst, method='average') ## average
## perform hierarchical clustering
dend <- color_branches(dend, k=3, col=c(lh.col, lc.col, se.col))
dend <- click_rotate(dend, plot=TRUE, plot_after=TRUE, continue=TRUE)
pdf('/Users/Avril/Desktop/dendro.pdf', width=4, height=4)
par(mar=c(2.1, 4.1, 2.1, 2.1))
dend %>%
  set('branches_lwd', 2.5) %>%
  set('labels_col', c(lc.col,lc.col,lc.col,lc.col,lc.col,lc.col,lc.col,lc.col,lc.col,se.col,se.col,se.col,se.col,se.col,lh.col,lh.col,lh.col,lh.col,lh.col)) %>%
  set('labels_cex', 0.5) %>%
  set('branches_k_color', k=3, value=c(lc.col, se.col, lh.col)) %>%
  # set('by_labels_branches_col', labels(dend)[1:14], type='all', TF_values=c('#956397',1)) %>% ## to make purple branch
  # raise.dendrogram(-22) %>%
  plot
dev.off()

##### Calculating AF values #####
## combine LC AFs calculated from individual genotypes and SE/LH AFs calculated from combining pool AFs
se.afs <- all.afs[,c(1:5)]
lh.afs <- all.afs[,c(6:10)]
## calculate population AFs from pool AFs for SE/LH
for(i in 1:nrow(se.afs)){
  se.afs$se.ref.af[i] <- ((sum(se.afs[i,c(1:5)])) / (length(which(!is.na(se.afs[i,c(1:5)])))))
}
for(i in 1:nrow(lh.afs)){
  lh.afs$lh.ref.af[i] <- ((sum(lh.afs[i,c(1:5)])) / (length(which(!is.na(lh.afs[i,c(1:5)])))))
}
## get LC REF AF calculated from individual genotypes
lc.ref.af <- lc.afs[which(rownames(lc.afs) %in% rownames(se.afs)),]
lc.ref.af <- lc.ref.af[order(rownames(lc.ref.af)),]
all(rownames(lc.ref.af) == rownames(se.afs))
all(rownames(lc.ref.af) == rownames(lh.afs))
## combine 3-pop REF AFs
pop.afs <- cbind(as.numeric(paste(lc.ref.af$lc.ref.af)), se.afs$se.ref.af, lh.afs$lh.ref.af)
pop.afs <- as.data.frame(pop.afs)
rownames(pop.afs) <- rownames(lc.ref.af)
colnames(pop.afs) <- c('lc.ref.af','se.ref.af','lh.ref.af')
pop.afs$chrom <- do.call(rbind, strsplit(rownames(pop.afs), split='|', fixed=TRUE))[,1]
pop.afs$pos <- do.call(rbind, strsplit(rownames(pop.afs), split='|', fixed=TRUE))[,2]
# write.table(pop.afs, '/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/snps_and_gene_expression/variant_calling/20_12_15_testing_scripts_for_ms_upload/pop_afs_23291_SNPs.txt', sep='\t', quote=FALSE)

##### Format population AF data for Fst calculation with Jacob's script #####
## Note: per Jacob's recommendation, all depths will be set equal to the number of individuals 
## per population, because for all loci, read depth exceeds the number of individuals sequenced;
## write file with necessary information for fst_from_pooled_freqs.pl
lc.dp <- seq(from=36, to=36, length.out=nrow(pop.afs))
pool.dp <- seq(from=20, to=20, length.out=nrow(pop.afs))
fst.input <- cbind(rownames(pop.afs), pop.afs$lc.ref.af, lc.dp, pop.afs$se.ref.af, pool.dp, pop.afs$lh.ref.af, pool.dp)
colnames(fst.input) <- c('Allele','Freq1','Size1','Freq2','Size2','Freq3','Size3')
write.table(fst.input, '/Users/Avril/Desktop/allele_freqs.txt', row.names=FALSE, sep='\t', quote=FALSE)
# >>> # run Jacob's script: perl fst_from_pooled_freqs.pl -f allele_freqs.txt -t -1
# <<< # read in output from that script
fsts.pool <- read.table('/Users/Avril/Desktop/ch3_working/High_Fsts_allele_freqs.txt', sep='\t', header=TRUE)
## write file of FST values to use for kNN analysis
colnames(fsts.pool) <- c('Allele','1,2','1,3','2,3')
write.table(fsts.pool, '/Users/Avril/Desktop/jacob_fsts.txt', quote=FALSE, sep='\t', row.names=FALSE)

colnames(fsts.pool) <- c('loc','fst.lc.se','fst.lc.lh','fst.se.lh')
pop.afs$loc <- rownames(pop.afs)
pop.afs <- merge(x=pop.afs, y=fsts.pool, by='loc')

##### z-transform Fst values and generate Manhattan plots #####
## calculate chromosome-specific z(Fst) values within chromosome
OUT <- NULL
for(i in unique(pop.afs$chrom)){
  sub <- pop.afs[which(pop.afs$chrom == i),]
  ## separate transformation
  sub$z.fsts.lc.se <- scale(sub$fst.lc.se, center=TRUE, scale=TRUE)
  sub$z.fsts.lc.lh <- scale(sub$fst.lc.lh, center=TRUE, scale=TRUE)
  sub$z.fsts.se.lh <- scale(sub$fst.se.lh, center=TRUE, scale=TRUE)
  
  OUT <- rbind(OUT, sub)
}
pop.afs <- OUT
rm(OUT)
#
##### Identify and annotate z(Fst) outliers #####
hist(pop.afs$z.fsts.lc.se, col=lc.vs.se, main='LC vs SE', xlab='z(FST)')
hist(pop.afs$z.fsts.lc.lh, col=lc.vs.lh, main='LC vs LH', xlab='z(FST)')
hist(pop.afs$z.fsts.se.lh, col=se.vs.lh, main='SE vs LH', xlab='z(FST)')

dim(pop.afs[which(pop.afs$z.fsts.lc.se >= 3),]) ## 501
dim(pop.afs[which(pop.afs$z.fsts.lc.lh >= 3),]) ## 366
dim(pop.afs[which(pop.afs$z.fsts.se.lh >= 3),]) ## 350
dim(pop.afs[which(pop.afs$z.fsts.lc.se >= 3 | 
              pop.afs$z.fsts.lc.lh >= 3 |
              pop.afs$z.fsts.se.lh >= 3),]) ## 1020

#### create exhaustive list of SNPs you might ever want annotations for:
## -- 231,915 LC SNPs
## -- 123,670 SE/LH SNPs
## >>> ## filter SnpEff output on cluster
## <<< ## read in filtered output (as a table, not as a vcf)
## read in filtered SnpEff annotations for LC SNPs
lc.snpeff <- read.table('/Users/Avril/Desktop/.scratch/analyses_ssalar_rna_seq/snpeff/LC_231915_SNPs/SPLIT_FIELDS_LC_231915_SNPs.ann.vcf',
                        sep='\t', header=TRUE)
se.lh.snpeff <- read.table('/Users/Avril/Desktop/.scratch/analyses_ssalar_rna_seq/snpeff/SE_LH_123670_SNPs/SPLIT_FIELDS_SE_LH_123670_SNPs.ann.vcf',
                           sep='\t', header=TRUE)

##### Compare kNN results to Fst and z(Fst) results -- generate illustrative plots #####
## read in weighted kNN scores
knnw.scores <- read.table('/Users/Avril/Desktop/all_ch_3_dirs/popgenome/pool_fsts_knnw_scores.txt', sep='\t', header=TRUE, check.names=FALSE)
rownames(knnw.scores) <- paste0(knnw.scores$chrom,'|',knnw.scores$pos)
## read in delta-Fst calculation results
## (calculated for 95% quantile outliers - can be pared down to any quantile cutoff)
delta.res <- read.table('/Users/Avril/Desktop/all_ch_3_dirs/popgenome/delta_Fst_results_95q.txt', sep='\t', check.names=FALSE, row.names=1)
nrow(delta.res)/nrow(knnw.scores) ## should be top 5% of kNN scores

## Manhattan plot highlighting kNN outliers
min.zfst.x <- min(c(pop.afs$z.fsts.lc.se, pop.afs$z.fsts.lc.lh, pop.afs$z.fsts.se.lh))
max.zfst.x <- max(c(pop.afs$z.fsts.lc.se, pop.afs$z.fsts.lc.lh, pop.afs$z.fsts.se.lh))

## highlight top 5% of kNN SNPs
# hilite <- rownames(delta.res)
## highlight top 1% of kNN SNPs
knnw.scores <- knnw.scores[order(-knnw.scores$knnw_scores),] ## sort by kNN score
n.q1 <- floor(nrow(knnw.scores)*0.01)
hilite <- rownames(knnw.scores)[1:n.q1]

## sort back for proper color assignments
knnw.scores <- knnw.scores[order(knnw.scores$chrom, knnw.scores$pos),]
pop.afs <- pop.afs[order(pop.afs$chrom, as.numeric(pop.afs$pos)),]
all(rownames(knnw.scores) == pop.afs$loc)
## prepare to shade outlier points by kNN score
ramp.pal <- colorRampPalette(c('yellow', 'red')) ## or lc.vs.se and magenta
# pop.afs$colour <- ramp.pal(30)[as.numeric(cut(knnw.scores$knnw_scores, breaks=30))] ## assign colors based on entire range of kNN scores
pop.afs[which(pop.afs$loc %in% hilite),'colour'] <- ramp.pal(30)[as.numeric(cut(knnw.scores[which(rownames(knnw.scores) %in% hilite),'knnw_scores'], breaks=30))] ## assign colors based on range of outlier kNN scores
## define colors, depending on whether MT chrom will be plotted
hi.cols <- pop.afs[which(pop.afs$loc %in% hilite),'colour']

all.chroms <- pop.afs
all.chroms$chrom <- as.factor(all.chroms$chrom)
all.chroms$chrom <- as.numeric(all.chroms$chrom)
all.chroms$pos <- as.numeric(all.chroms$pos)
lowlim <- min.zfst.x
uplim <- max.zfst.x

##### Manhattan plotting - LC vs SE #####
colour <- '#47304a'
filler <- 'gray10'
comp <- 'z.fsts.lc.se'

## can turn on for unique ylim for each plot
lowlim <- min(all.chroms[,comp])
uplim <- max(all.chroms[,comp])

## manhattan plot using ggrepel package
don <- all.chroms %>%
  group_by(chrom) %>%
  summarise(chr_len=max(pos)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(all.chroms, ., by=c('chrom'='chrom')) %>%
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot) %>%
  # mutate(is_highlight=ifelse(loc %in% filt.outliers$X, 'yes', 'no')) %>%
  mutate(is_highlight=ifelse(loc %in% hilite, 'yes', 'no')) %>%
  mutate(is_annotate=ifelse(loc %in% hilite, 'yes', 'no'))

hi.cols <- subset(don, is_highlight=='yes' & don[,comp] >= 3)[,'colour']

axisdf <- don %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

pdf(paste0('/Users/Avril/Desktop/zFst_SEP_',comp,'.pdf'), width=10, height=2.5)
# Make the plot
ggplot(don, aes(x=BPcum, y=z.fsts.lc.se)) +
  # Show all points
  geom_point( aes(color=as.factor(chrom)), alpha=0.9, size=1.3) + #### CAN CHANGE ALPHA WHEN HIGHLIGHTING OUTLIER SNPS
  scale_color_manual(values = rep(c(filler, colour), 22 )) +
  # custom X axis:
  scale_x_continuous( expand = c(0.01,0), label = c(axisdf$chrom[c(1,5,10,15,20,25)],'scaffolds'), breaks= axisdf$center[c(1,5,10,15,20,25,31)],
                      name = 'Chromosome') + ## to label every 5 chrom + MT + scafs
  scale_y_continuous(expand = c(0, 0) , limits=c(-0.5, 8.5), name = bquote(z(italic(F[ST]))), breaks=c(0,2.5,5,7.5)) +
  ## add z(Fst) cutoff line
  geom_hline(yintercept=3, linetype='solid', color='red', size=0.75) +
  
  ## ------------------------- to highlight kNN outlier SNPs
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes" & don[,comp] >= 3), color='black', size=2.5, pch=23,
             bg=hi.cols) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=hlite.names), size=5) +
  ## -------------------------

# Custom the theme:
theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.y = element_line(colour = 'black', size=0.25, linetype='solid'),
    axis.line.x = element_line(colour = 'black', size=0.25, linetype='solid'),
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)) ## LC/SE = 5
  )
dev.off()

##### Manhattan plotting - LC vs LH #####
colour <- '#229284'
filler <- 'gray10'
comp <- 'z.fsts.lc.lh'

## can turn on for unique ylim for each plot
lowlim <- min(all.chroms[,comp])
uplim <- max(all.chroms[,comp])

## manhattan plot using ggrepel package
don <- all.chroms %>%
  group_by(chrom) %>%
  summarise(chr_len=max(pos)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(all.chroms, ., by=c('chrom'='chrom')) %>%
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot) %>%
  # mutate(is_highlight=ifelse(loc %in% filt.outliers$X, 'yes', 'no')) %>%
  mutate(is_highlight=ifelse(loc %in% hilite, 'yes', 'no')) %>%
  mutate(is_annotate=ifelse(loc %in% hilite, 'yes', 'no'))

hi.cols <- subset(don, is_highlight=='yes' & don[,comp] >= 3)[,'colour']

axisdf <- don %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# pdf(paste0('/Users/Avril/Desktop/zFst_SEP_',comp,'.pdf'), width=12, height=4)
pdf(paste0('/Users/Avril/Desktop/zFst_SEP_',comp,'.pdf'), width=10, height=2.5)
# Make the plot
ggplot(don, aes(x=BPcum, y=z.fsts.lc.lh)) +
  # ggplot(don, aes(x=BPcum, y=fst.se.lh)) +  ## try plotting with just Fst instead of z(Fst)
  # Show all points
  geom_point( aes(color=as.factor(chrom)), alpha=0.9, size=1.3) + #### CAN CHANGE ALPHA WHEN HIGHLIGHTING OUTLIER SNPS
  scale_color_manual(values = rep(c(filler, colour), 22 )) +
  # custom X axis:
  # scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_x_continuous( expand = c(0.01,0), label = c(axisdf$chrom[c(1,5,10,15,20,25)],'scaffolds'), breaks= axisdf$center[c(1,5,10,15,20,25,31)],
                      name = 'Chromosome') + ## to label every 5 chrom + MT + scafs
  scale_y_continuous(expand = c(0, 0) , limits=c(-0.5, 6), name = bquote(z(italic(F[ST]))), breaks=c(0,2,4,6)) +
  ## add z(Fst) cutoff line
  geom_hline(yintercept=3, linetype='solid', color='red', size=0.75) +
  
  ## ------------------------- to highlight kNN outlier SNPs
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes" & don[,comp] >= 3), color='black', size=2.5, pch=23,
             bg=hi.cols) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=hlite.names), size=5) +
  ## -------------------------

# Custom the theme:
theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.y = element_line(colour = 'black', size=0.25, linetype='solid'),
    axis.line.x = element_line(colour = 'black', size=0.25, linetype='solid'),
    axis.title.y = element_text(margin = margin(t = 0, r = 12.4, b = 0, l = 0)) 
  )
dev.off()

##### Manhattan plotting - SE vs LH #####
colour <- '#e27b0f'
filler <- 'gray10'
comp <- 'z.fsts.se.lh'

## can turn on for unique ylim for each plot
lowlim <- min(all.chroms[,comp])
uplim <- max(all.chroms[,comp])

## manhattan plot using ggrepel package
don <- all.chroms %>%
  group_by(chrom) %>%
  summarise(chr_len=max(pos)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(all.chroms, ., by=c('chrom'='chrom')) %>%
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot) %>%
  # mutate(is_highlight=ifelse(loc %in% filt.outliers$X, 'yes', 'no')) %>%
  mutate(is_highlight=ifelse(loc %in% hilite, 'yes', 'no')) %>%
  mutate(is_annotate=ifelse(loc %in% hilite, 'yes', 'no'))

hi.cols <- subset(don, is_highlight=='yes' & don[,comp] >= 3)[,'colour']

axisdf <- don %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# pdf(paste0('/Users/Avril/Desktop/zFst_SEP_',comp,'.pdf'), width=12, height=4)
pdf(paste0('/Users/Avril/Desktop/zFst_SEP_',comp,'.pdf'), width=10, height=2.5)
# Make the plot
ggplot(don, aes(x=BPcum, y=z.fsts.se.lh)) +
  # ggplot(don, aes(x=BPcum, y=fst.se.lh)) +  ## try plotting with just Fst instead of z(Fst)
  # Show all points
  geom_point( aes(color=as.factor(chrom)), alpha=0.9, size=1.3) + #### CAN CHANGE ALPHA WHEN HIGHLIGHTING OUTLIER SNPS
  scale_color_manual(values = rep(c(filler, colour), 22 )) +
  # custom X axis:
  # scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_x_continuous( expand = c(0.01,0), label = c(axisdf$chrom[c(1,5,10,15,20,25)],'scaffolds'), breaks= axisdf$center[c(1,5,10,15,20,25,31)],
                      name = 'Chromosome') + ## to label every 5 chrom + MT + scafs
  scale_y_continuous(expand = c(0, 0) , limits=c(-0.5, 6), name = bquote(z(italic(F[ST]))), breaks=c(0,2,4,6)) +
  ## add z(Fst) cutoff line
  geom_hline(yintercept=3, linetype='solid', color='red', size=0.75) +
  
  ## ------------------------- to highlight kNN outlier SNPs
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes" & don[,comp] >= 3), color='black', size=2.5, pch=23,
             bg=hi.cols) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=hlite.names), size=5) +
  ## -------------------------

# Custom the theme:
theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.y = element_line(colour = 'black', size=0.25, linetype='solid'),
    axis.line.x = element_line(colour = 'black', size=0.25, linetype='solid'),
    axis.title.y = element_text(margin = margin(t = 0, r = 12.4, b = 0, l = 0)) 
  )
dev.off()

##### Annotating kNN outliers #####
colnames(lc.snpeff) <- c('chrom','pos','ref','alt','allele','effect','impact','gene.name','gene.ID','feature.type','feature.ID',
                         'biotype','rank.or.total','hgvs.c','hgvs.p','cdna.pos','len.cdna','pos.aa','len.aa','dist.to.feat','errors')
## read in reference annotation for matching up gene names
full.assembly <- readGFF('/Users/Avril/Desktop/all_ch_3_dirs/02_output_tables_all_analyses/GCF_000233375.1_ICSASG_v2_genomic.gff')
## read in DEG + AP lists
treat.degs <- read.csv('/Users/Avril/Desktop/all_ch_3_dirs/02_output_tables_all_analyses/3616_degs.csv')
fam.degs <- read.csv('/Users/Avril/Desktop/all_ch_3_dirs/02_output_tables_all_analyses/1446_genes.csv')
ap.genes <- read.table('/Users/Avril/Desktop/all_ch_3_dirs/02_output_tables_all_analyses/unique_a_priori_gene_list.txt')
## read in gene description information
descrip <- read.csv('/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/stringtie_featureCounts_deseq2_DGE_analyses/ICSASG_v2_gene_list.csv')
## subset reference annotation to link up SnpEff output column 'gene.ID' and ref annotation column 'ID'
sub.ass <- as.data.frame(full.assembly[,c(1:5,9,10,21,22,23,25,26)])

## filtering SnpEff annotations
## limit by kNN outliers first
## highlight top 1% of kNN SNPs
knnw.scores <- knnw.scores[order(-knnw.scores$knnw_scores),] ## sort by kNN score
n.q1 <- floor(nrow(knnw.scores)*0.01) ## can set quantile cutoff here
hilite <- rownames(knnw.scores)[1:n.q1] ## n=232 SNPs if quantile = top 1%

lc.snpeff$loc <- paste0(lc.snpeff$chrom,'|',lc.snpeff$pos)
lc.snpeff <- lc.snpeff[which(lc.snpeff$loc %in% hilite),]
length(unique(lc.snpeff$loc)) ## 225 SNPs -- 7 w/ no annotations

## limit to just  zFst outliers
knn.out <- pop.afs[which(pop.afs$loc %in% hilite),]
knn.out <- knn.out[which(knn.out$z.fsts.lc.se >= 3 | knn.out$z.fsts.lc.lh >= 3 | knn.out$z.fsts.se.lh >= 3),] ## 209 kNN outliers with z(FST) >=3
## quantify outliers for each comparison
nrow(knn.out[which(knn.out$z.fsts.lc.se >= 3),]) ## 37
nrow(knn.out[which(knn.out$z.fsts.lc.lh >= 3),]) ## 160
nrow(knn.out[which(knn.out$z.fsts.se.lh >= 3),]) ## 145
## limit to annotated SNPs
knn.out <- knn.out[which(knn.out$loc %in% lc.snpeff$loc),]
knn.out <- knn.out[,c(1,10:12)]
dim(knn.out) ## 203 SNPs

lc.snpeff <- lc.snpeff[which(lc.snpeff$loc %in% knn.out$loc),]
length(unique(lc.snpeff$loc)) ## 203 SNPs
## feature types for these snps?
table(lc.snpeff$feature.type) ## largely transcripts, but also 27 intergenic_region;
## first focus on SNPs in transcripts (then on SNPs in intergenic regions)
## NARROW IT DOWN TO SNPs IN TRANSCRIPTS FIRST
trans.snps <- lc.snpeff[which(lc.snpeff$feature.type=='transcript' & lc.snpeff$biotype != 'pseudogene'),] ## added in pseudogene requirement on 9/9/20,
## because I found that some SNPs were annotated to pseudogenes AND protein-coding genes, but only the pseudogene annotations
## were retained in the final annotation set
length(unique(trans.snps$loc)) ## 196 / 209 (94%) SNPs are covered in this category

## remove all annotations with an error associated
trans.snps <- trans.snps[which(trans.snps$errors==''),] 
length(unique(trans.snps$loc)) ## --> 190 / 209 SNPs

## merge on sub.ass$transcript_id and trans.snps$feature.ID
ids.genes <- sub.ass[,c('transcript_id', 'gene')]
colnames(ids.genes) <- c('feature.ID','gene')
ids.genes <- ids.genes[!duplicated(ids.genes),]
## add parent gene information to transcript SNPs using ID:gene key
dim(merge(x=trans.snps, y=ids.genes, by='feature.ID', all.x=TRUE))
trans.snps <- merge(x=trans.snps, y=ids.genes, by='feature.ID', all.x=TRUE)
trans.snps <- trans.snps[order(trans.snps$chrom, trans.snps$pos),]

## save annotations with gene labels and delete from main dataframe
table(is.na(trans.snps$gene)) ## this only labeled 0 / 610 annotations
OUT <- NULL
OUT <- rbind(OUT, trans.snps[!is.na(trans.snps$gene),]) ## save annotations in OUT
trans.snps <- trans.snps[,-23] ## get rid of empty gene column

## try merging using sub.ass$Parent to label exons/transcript annotations
ids.genes <- sub.ass[,c('Parent', 'gene')]
colnames(ids.genes) <- c('feature.ID','gene')
ids.genes <- ids.genes[!duplicated(ids.genes),]
## add parent gene information to transcript SNPs using ID:gene key
dim(merge(x=trans.snps, y=ids.genes, by='feature.ID', all.x=TRUE))
trans.snps <- merge(x=trans.snps, y=ids.genes, by='feature.ID', all.x=TRUE)
trans.snps <- trans.snps[order(trans.snps$chrom, trans.snps$pos),]

## save annotations with gene labels and delete from main dataframe
table(is.na(trans.snps$gene)) ## 0 genes unassigned
OUT <- rbind(OUT, trans.snps[!is.na(trans.snps$gene),]) ## add annotations to OUT
trans.snps <- trans.snps[,-23] ## get rid of gene column

## which loc still don't have a gene associated?
trans.snps[which(trans.snps$loc %notin% OUT$loc),] ## --> 0 SNPs in trans.snps unannotated
trans.snps <- OUT ## replace trans.snps with completed annotations
length(unique(trans.snps$loc)) ## 190 SNPs have annotations without annotations

length(unique(lc.snpeff[which(lc.snpeff$loc %notin% trans.snps$loc),'loc'])) ## 13 SNPs
## examine unannotated SNPs
unannot <- lc.snpeff[which(lc.snpeff$loc %notin% trans.snps$loc),]
unannot <- unannot[which(unannot$errors==''),]
table(unannot$feature.type) ## these SNPs either occur in intergenic regions or have no error-free annotations; 
unannot[unannot$feature.type == 'transcript',] ## 'transcripts' are annotated to pseudogenes

## loop over each SNP, only keeping 1 annotation per
## SNP:gene association, unless there are multiple effects for the same assocation (i.e., transcript-
## specific effects)
OUT <- NULL
OUT1 <- NULL
for(i in 1:length(unique(trans.snps$loc))){
  temp <- trans.snps[which(trans.snps$loc == unique(trans.snps$loc)[i]),]
  if(length(unique(temp$gene)) == 1 & ## if there is only 1 gene associated with this loc &
     length(unique(temp$effect)) == 1){ ## only 1 effect observed for this association,
    OUT <- rbind(OUT, temp[1,]) ## save that assocation + effect combination
  }
  else{ ## otherwise, 
    for(j in unique(temp$gene)){ ## loop over each gene associated with the loc,
      sub.temp <- temp[which(temp$gene == j),] ## and subset by gene.
      if(length(unique(sub.temp$effect)) == 1){ ## if there is only one effect for that gene,
        OUT <- rbind(OUT, sub.temp[1,]) ## save it.
      }
      else{
        OUT1 <- rbind(OUT1, sub.temp) ## otherwise, save SNP:gene assoc. for further filtering below
      }
    }
  }
  print(i/length(unique(trans.snps$loc)))
}

## for each SNP in OUT1, keep 1 instance of each type of effect
OUT2 <- NULL
for(i in 1:length(unique(OUT1$loc))){
  temp <- OUT1[which(OUT1$loc == unique(OUT1$loc)[i]),]
  temp <- temp[!duplicated(temp$effect),]
  OUT2 <- rbind(OUT2, temp)
  print(i/length(unique(OUT1$loc)))
}

## I think this is as narrowed-down as it can get, as far as in-/proximal-to-transcript
## SNPs go;
## OUT = unique SNP:gene associations; if a SNP is associated with >1 transcript but the effects are the same for all,
## just the first annotation is saved per SNP:gene association.
## OUT2 = for each SNP:gene association with >1 effect type, keep only 1 instance of each effect type.
## Combine these two datasets:
table(unique(trans.snps$loc) %in% unique(c(OUT$loc, OUT2$loc))) ## 196 SNPs w/ error-free annotations
trans.snps <- rbind(OUT, OUT2)
trans.snps <- trans.snps[order(trans.snps$chrom, trans.snps$pos),] ## final, filtered list of annotations for SNPs in/near transcripts
## add gene description information
sub.desc <- descrip[,c(6,8)]
sub.desc <- sub.desc[!duplicated(sub.desc$Symbol),]
colnames(sub.desc)[1] <- 'gene'
trans.snps <- merge(x=trans.snps, y=sub.desc, by='gene', all.x=TRUE)
# ## add delta-Fst values for gene interpretation
# temp <- delta.res[,c(1:3,5)]
# trans.snps <- merge(x=trans.snps, y=temp, by='loc', all.x=TRUE)
## add in DEG information (no AP genes in kNN 1% outliers)
trans.snps[which(trans.snps$gene %in% fam.degs$gene.name),'FAM.DEG'] <- TRUE
trans.snps[which(trans.snps$gene %in% treat.degs$gene.name),'TREAT.DEG'] <- TRUE

write.table(trans.snps, '/Users/Avril/Desktop/kNN_190_outlier_SNPs_q_01_annotations.txt', sep='\t', quote=FALSE, row.names=FALSE)
## write gene list
write.table(unique(trans.snps$gene), '/Users/Avril/Desktop/kNN_outliers.txt', quote=FALSE, sep='\t', row.names=FALSE)

knn.out <- pop.afs[which(pop.afs$loc %in% trans.snps$loc),]
knn.out <- knn.out[,c(1,10:12)]
colnames(knn.out)[2:4] <- c('z.fsts.lc.se','z.fsts.lc.lh','z.fsts.se.lh')
knn.out$zfst.out <- knn.out$z.fsts.lc.se >= 3 | knn.out$z.fsts.lc.lh >= 3 | knn.out$z.fsts.se.lh >= 3

trans.snps <- merge(x=trans.snps, y=knn.out, by='loc', fixed=TRUE, all.x=TRUE)
# write.table(trans.snps, '/Users/Avril/Desktop/kNN_252_outlier_SNPs_q_01_annotations_w_zFst_info.txt', sep='\t', quote=FALSE, row.names=FALSE)
trans.snps <- trans.snps[which(trans.snps$zfst.out == TRUE),]
# dim(trans.snps)
# length(unique(trans.snps$loc))
## add comparison-specific information
trans.snps[trans.snps$z.fsts.lc.lh >= 3 & trans.snps$z.fsts.se.lh >= 3, 'anad'] <- TRUE
trans.snps[trans.snps$z.fsts.lc.lh >= 3 & trans.snps$z.fsts.lc.se >= 3, 'hatch'] <- TRUE
write.table(trans.snps, '/Users/Avril/Desktop/kNN_190_outlier_SNPs_q_01_annotations_w_zFst_info_PSEUDOGENE_UPDATE.txt', sep='\t', quote=FALSE, row.names=FALSE)

##### Calculating pooled heterozygosity #####
## SE
se.depths <- pools.10.dp[,c(25,3:9,15:19)]
se.depths$tot.ref <- rowSums(se.depths[,c(4:8)])
se.depths$tot.alt <- rowSums(se.depths[,c(9:13)])-se.depths$tot.ref
se.depths$chrom <- do.call(rbind, strsplit(as.character(se.depths$loc), split='|', fixed=TRUE))[,1]
se.depths <- se.depths[,c(1:3,14,15,16)]
for(i in 1:nrow(se.depths)){
  if(se.depths$tot.ref[i] > se.depths$tot.alt[i]){
    se.depths$maj[i] <- se.depths$tot.ref[i]
    se.depths$min[i] <- se.depths$tot.alt[i]
  }
  else{
    se.depths$maj[i] <- se.depths$tot.alt[i]
    se.depths$min[i] <- se.depths$tot.ref[i]
  }
}

## for each SNP
se.depths$se.per.loc <- (2 * (se.depths$maj) * (se.depths$min)) / ((se.depths$maj) + (se.depths$min))^2
se.hp <- mean(se.depths$se.per.loc)

## LH
lh.depths <- pools.10.dp[,c(25,10:14,20:24)]
lh.depths$tot.ref <- rowSums(lh.depths[,c(2:6)])
lh.depths$tot.alt <- rowSums(lh.depths[,c(7:11)])-lh.depths$tot.ref
lh.depths$chrom <- do.call(rbind, strsplit(as.character(lh.depths$loc), split='|', fixed=TRUE))[,1]
lh.depths <- lh.depths[,c(1,12,13,14)]
for(i in 1:nrow(lh.depths)){
  if(lh.depths$tot.ref[i] > lh.depths$tot.alt[i]){
    lh.depths$maj[i] <- lh.depths$tot.ref[i]
    lh.depths$min[i] <- lh.depths$tot.alt[i]
  }
  else{
    lh.depths$maj[i] <- lh.depths$tot.alt[i]
    lh.depths$min[i] <- lh.depths$tot.ref[i]
  }
}

## for each SNP
lh.depths$lh.per.loc <- (2 * (lh.depths$maj) * (lh.depths$min)) / ((lh.depths$maj) + (lh.depths$min))^2
lh.hp <- mean(lh.depths$lh.per.loc)

## LC, treating all 36 samples as a pool by summing allele read depth across them all
all(rownames(total.dp) == rownames(ref.dp))
lc.depths <- as.data.frame(cbind(rownames(total.dp), rowSums(ref.dp), (rowSums(total.dp)-rowSums(ref.dp))))
colnames(lc.depths) <- c('loc','tot.ref','tot.alt')
lc.depths <- lc.depths[which(lc.depths$loc %in% pop.afs$loc),]
lc.depths$tot.ref <- as.numeric(paste(lc.depths$tot.ref))
lc.depths$tot.alt <- as.numeric(paste(lc.depths$tot.alt))
lc.depths$chrom <- do.call(rbind, strsplit(as.character(lc.depths$loc), split='|', fixed=TRUE))[,1]
for(i in 1:nrow(lc.depths)){
  if(lc.depths$tot.ref[i] > lc.depths$tot.alt[i]){
    lc.depths$maj[i] <- lc.depths$tot.ref[i]
    lc.depths$min[i] <- lc.depths$tot.alt[i]
  }
  else{
    lc.depths$maj[i] <- lc.depths$tot.alt[i]
    lc.depths$min[i] <- lc.depths$tot.ref[i]
  }
}

## for each SNP
lc.depths$lc.per.loc <- (2 * (lc.depths$maj) * (lc.depths$min)) / ((lc.depths$maj) + (lc.depths$min))^2
lc.hp <- mean(lc.depths$lc.per.loc)

# pdf('/Users/Avril/Desktop/hp_whole_dataset.pdf', width=3, height=4)
barplot(c(lc.hp, se.hp, lh.hp), axis.lty=1, xpd=FALSE,
        names.arg=c('LC','SE','LH'), col=c(lc.col, se.col, lh.col),
        ylab=bquote(italic(H[P])), cex.axis=0.8, cex.names=0.8)
# dev.off()

## calculate per chromosome?
## SE
se.chrom.hp <- NULL
for(i in unique(se.depths$chrom)){
  sub <- se.depths[which(se.depths$chrom == i),]
  # temp <- c(((2 * sum(sub$maj) * sum(sub$min)) / (sum(sub$maj) + sum(sub$min))^2), i)
  temp <- mean(sub$se.per.loc)
  se.chrom.hp <- rbind(se.chrom.hp, temp)
}

## LH
lh.chrom.hp <- NULL
for(i in unique(lh.depths$chrom)){
  sub <- lh.depths[which(lh.depths$chrom == i),]
  # temp <- c(((2 * sum(sub$maj) * sum(sub$min)) / (sum(sub$maj) + sum(sub$min))^2), i)
  temp <- mean(sub$lh.per.loc)
  lh.chrom.hp <- rbind(lh.chrom.hp, temp)
}

## LC
lc.chrom.hp <- NULL
for(i in unique(lc.depths$chrom)){
  sub <- lc.depths[which(lc.depths$chrom == i),]
  # temp <- c(((2 * sum(sub$maj) * sum(sub$min)) / (sum(sub$maj) + sum(sub$min))^2), i)
  temp <- mean(sub$lc.per.loc)
  lc.chrom.hp <- rbind(lc.chrom.hp, temp)
}

## plot all together
x.lc.chrom.hp <- cbind(seq(from=1, to=1, length=nrow(lc.chrom.hp)), lc.chrom.hp, lc.col)
x.se.chrom.hp <- cbind(seq(from=2, to=2, length=nrow(se.chrom.hp)), se.chrom.hp, se.col)
x.lh.chrom.hp <- cbind(seq(from=3, to=3, length=nrow(lh.chrom.hp)), lh.chrom.hp, lh.col)
all.chrom.hp <- rbind(x.lc.chrom.hp, x.se.chrom.hp, x.lh.chrom.hp)
# pdf('/Users/Avril/Desktop/hp_by_chrom.pdf', width=3, height=5)
par(mgp=c(2.3,1,0))
plot(jitter(as.numeric(all.chrom.hp[,1]), factor=0.7), all.chrom.hp[,2], pch=19, col=alpha(c(all.chrom.hp[,3]), 0.6),
     ylab=bquote(italic(H[P])), xaxt='n', ylim=c(0.175, 0.35), xlim=c(0.5, 3.5), xlab='Population', yaxt='n')
  lines(x=c(0.6, 1.4), y=c(lc.hp, lc.hp), col=lc.col, lwd=4)
  lines(x=c(1.6, 2.4), y=c(se.hp, se.hp), col=se.col, lwd=4)
  lines(x=c(2.6, 3.4), y=c(lh.hp, lh.hp), col=lh.col, lwd=4)
# ## optional: add means calc from chrom values
# lines(x=c(0.6, 1.4), y=c(mean(as.numeric(lc.chrom.hp[,1])), mean(as.numeric(lc.chrom.hp[,1]))), col=lc.col, lwd=4, lty=3)
# lines(x=c(1.6, 2.4), y=c(mean(as.numeric(se.chrom.hp[,1])), mean(as.numeric(se.chrom.hp[,1]))), col=se.col, lwd=4, lty=3)
# lines(x=c(2.6, 3.4), y=c(mean(as.numeric(lh.chrom.hp[,1])), mean(as.numeric(lh.chrom.hp[,1]))), col=lh.col, lwd=4, lty=3)
  axis(1, at=c(1,2,3), labels=c('LC','SE','LH'))
  axis(2, at=c(0.20,0.25,0.30,0.35), labels=c('0.20','0.25','0.30','0.35'))
# dev.off()
  
##### Randomization test: differences among populations in Hp #####
lc.pool <- as.numeric(lc.chrom.hp[,1])
se.pool <- as.numeric(se.chrom.hp[,1])
lh.pool <- as.numeric(lh.chrom.hp[,1])

pdf('/Users/Avril/Desktop/hp_rando_tests.pdf', width=6, height=12)
par(mfrow=c(3,1), mar=c(5.1,5.1,4.1,2.1))
OUT <- NULL
for(i in 1:10000){
  temp <- sample(c(lc.pool, se.pool), size=29)
  OUT <- c(OUT, mean(temp))
  print(i)
}
hist(OUT, col='grey', xlim=c(min(c(mean(lc.pool), mean(se.pool)))-0.005, max(c(mean(lc.pool), mean(se.pool)))+0.005),
     ylim=c(0,2000), main='LC vs SE', xlab=expression('Mean sample '*italic(H)[P]), cex.axis=1.5, cex.lab=1.5)
abline(v=mean(lc.pool), col=lc.col, lwd=2, lty=2)
abline(v=mean(se.pool), col=se.col, lwd=2, lty=2)

OUT <- NULL
for(i in 1:10000){
  temp <- sample(c(lc.pool, lh.pool), size=29)
  OUT <- c(OUT, mean(temp))
  print(i)
}
hist(OUT, col='grey', xlim=c(min(c(mean(lc.pool), mean(lh.pool)))-0.005, max(c(mean(lc.pool), mean(lh.pool)))+0.005),
     ylim=c(0,2500), main='LC vs LH', xlab=expression('Mean sample '*italic(H)[P]), cex.axis=1.5, cex.lab=1.5)
abline(v=mean(lc.pool), col=lc.col, lwd=2, lty=2)
abline(v=mean(lh.pool), col=lh.col, lwd=2, lty=2)

OUT <- NULL
for(i in 1:10000){
  temp <- sample(c(se.pool, lh.pool), size=29)
  OUT <- c(OUT, mean(temp))
  print(i)
}
# hist(OUT, col='grey', xlim=c(min(c(mean(se.pool), mean(lh.pool)))-0.01, max(c(mean(se.pool), mean(lh.pool)))+0.01), main='SE vs LH')
hist(OUT, col='grey', main='SE vs LH', ylim=c(0,1000), xlab=expression('Mean sample '*italic(H)[P]), cex.axis=1.5, cex.lab=1.5, breaks=30)
abline(v=mean(se.pool), col=se.col, lwd=2, lty=2)
abline(v=mean(lh.pool), col=lh.col, lwd=2, lty=2)
dev.off()

length(OUT[which(OUT > mean(se.pool))])/length(OUT) + length(OUT[which(OUT < mean(lh.pool))])/length(OUT) ## 0.057

##### Estimating effect of full-siblings on Hp in LC #####
## reorder LC genotype information to make sampling 1 indiv/family easier
total.dp <- total.dp[,order(colnames(total.dp))]
total.dp <- total.dp[,c(1:2,11:12,3:10,13:36)]
ref.dp <- ref.dp[,order(colnames(ref.dp))]
ref.dp <- ref.dp[,c(1:2,11:12,3:10,13:36)]

OUT.SAMPS <- NULL
OUT <- NULL
j <- 1
n <- 1001 ## how many iterations you want to run + 1
while(j < n){
  temp.total.dp <- total.dp
  temp.ref.dp <- ref.dp
  ## randomly subsample to get 1 indiv per family
  f.1 <- sample(1:4, size=1)
  f.11 <- sample(5:8, size=1)
  f.13 <- sample(9:12, size=1)
  f.3 <- sample(13:16, size=1)
  f.4 <- sample(17:20, size=1)
  f.6 <- sample(21:24, size=1)
  f.7 <- sample(25:28, size=1)
  f.8 <- sample(29:32, size=1)
  f.9 <- sample(33:36, size=1)
  
  temp.total.dp <- temp.total.dp[,c(f.1, f.11, f.13, f.3, f.4, f.6, f.7, f.8, f.9)]
  temp.ref.dp <- temp.ref.dp[,c(f.1, f.11, f.13, f.3, f.4, f.6, f.7, f.8, f.9)]
  OUT.SAMPS <- rbind(OUT.SAMPS, colnames(temp.total.dp))
  
  ## Calculating pooled heterozygosity ##
  lc.depths <- as.data.frame(cbind(rownames(temp.total.dp), rowSums(temp.ref.dp), (rowSums(temp.total.dp)-rowSums(temp.ref.dp))))
  colnames(lc.depths) <- c('loc','tot.ref','tot.alt')
  lc.depths$tot.ref <- as.numeric(paste(lc.depths$tot.ref))
  lc.depths$tot.alt <- as.numeric(paste(lc.depths$tot.alt))
  lc.depths$chrom <- do.call(rbind, strsplit(as.character(lc.depths$loc), split='|', fixed=TRUE))[,1]
  for(i in 1:nrow(lc.depths)){
    if(lc.depths$tot.ref[i] > lc.depths$tot.alt[i]){
      lc.depths$maj[i] <- lc.depths$tot.ref[i]
      lc.depths$min[i] <- lc.depths$tot.alt[i]
    }
    else{
      lc.depths$maj[i] <- lc.depths$tot.alt[i]
      lc.depths$min[i] <- lc.depths$tot.ref[i]
    }
  }
  lc.depths$lc.per.loc <- (2 * (lc.depths$maj) * (lc.depths$min)) / ((lc.depths$maj) + (lc.depths$min))^2
  loop.lc.hp <- mean(lc.depths$lc.per.loc)
  OUT <- c(OUT, loop.lc.hp)
  
  print(paste0(j,' - ',Sys.time()))
  j <- j+1
}
safe <- OUT
### save output of iteration in case need to run additional iterations/edit figure
# write.table(OUT, '/Users/Avril/Desktop/hp_LC_nonsibling_sampling.txt', sep='\t', quote=FALSE, row.names=FALSE)
OUT <- read.table('/Users/Avril/Desktop/ch3_working/hp_LC_nonsibling_sampling.txt', sep='\t', header=TRUE)
mean(OUT$x)
sd(OUT$x)
# pdf('/Users/Avril/Desktop/hp_whole_dataset.pdf', width=4, height=5)
par(mgp=c(3,1.5,0))
bp <- barplot(c(mean(OUT$x), lc.hp, se.hp, lh.hp), axis.lty=1, xpd=FALSE, ylim=c(0.20, 0.28),
              names.arg=c('LC sibling\nsampling','LC','SE','LH'), col=c(alpha(lc.col, 0.7), lc.col, se.col, lh.col),
              ylab=bquote(italic(H)[P]), cex.axis=0.8, cex.names=0.8, border=TRUE, xaxt='n')
  arrows(x0=bp[1], y0=mean(OUT$x)-sd(OUT$x), x1=bp[1], y1=mean(OUT$x)+sd(OUT$x), angle=90, code=3, length=0.075)
  axis(1, at=bp, labels=c('LC non-sibling\nsampling','LC','SE','LH'), lwd=0, lwd.ticks=1, cex.axis=0.8)
  abline(h=0.20, lwd=2)
# dev.off()

##### Estimating effect of full-siblings on FST in LC comparisons #####
## run scripts on cluster (10,000 iterations), read output in here to summarize
## read in results from cluster - 10,000 iterations
SE.FST.RANGE <- read.table('/Users/Avril/Desktop/.scratch/analyses_ssalar_rna_seq/lc_sibship_iterations/fst/summary_stats/10000_n9_lc_se_fst_mean_SD.txt', sep='\t', header=TRUE)
colnames(SE.FST.RANGE) <- c('mean','sd','min','max')
LH.FST.RANGE <- read.table('/Users/Avril/Desktop/.scratch/analyses_ssalar_rna_seq/lc_sibship_iterations/fst/summary_stats/10000_n9_lc_lh_fst_mean_SD.txt', sep='\t', header=TRUE)
colnames(LH.FST.RANGE) <- c('mean','sd','min','max')

# pdf('/Users/Avril/Desktop/n9_n36_hists.pdf', width=5, height=5)
hist(SE.FST.RANGE$mean, col=lc.vs.se, main='LC vs SE', xlab=expression('Mean '*italic(F)[ST]*' difference'), ylim=c(0,14000), xlim=c(0,0.20))
  abline(v=median(SE.FST.RANGE$mean), col='red', lty=2, lwd=1.5)
hist(LH.FST.RANGE$mean, col=lc.vs.lh, main='LC vs LH', xlab=expression('Mean '*italic(F)[ST]*' difference'), xlim=c(0,0.20))
  abline(v=median(LH.FST.RANGE$mean), col='red', lty=2, lwd=1.5)
# dev.off()

##### Write SI Tables #####
## S1: Outlier SNPs
## Columns: chrom, pos, knnw.score, z.fst.lc.se, z.fst.lc.lh, z.fst.se.lh
knnw.scores$loc <- rownames(knnw.scores)
s1 <- merge(x=pop.afs, y=knnw.scores, by='loc')
s1 <- s1[which(s1$loc %in% hilite),]
s1 <- s1[which(s1$z.fsts.lc.se >= 3 | s1$z.fsts.lc.lh >= 3 | s1$z.fsts.se.lh >= 3),]
s1 <- s1[,c('chrom.x','pos.x','knnw_scores','z.fsts.lc.se','z.fsts.lc.lh','z.fsts.se.lh')]
colnames(s1) <- c('chrom','pos','knnw.score','z.fst.lc.se','z.fst.lc.lh','z.fst.se.lh')
s1 <- s1[order(s1$chrom, as.numeric(s1$pos)),]
write.csv(s1, '/Users/Avril/Desktop/S1.csv', row.names=FALSE, quote=FALSE)

## S2: All SnpEff annotations for outlier SNPs
## Columns: chrom, pos, ref, alt, gene, description, effect, impact, z.fst.lc.se, z.fst.lc.lh, z.fst.se.lh
s2 <- trans.snps[which(trans.snps$loc %in% hilite),]
s2 <- s2[which(s2$z.fsts.lc.se >= 3 | s2$z.fsts.lc.lh >= 3 | s2$z.fsts.se.lh >= 3),]
s2 <- s2[,c('chrom','pos','ref','alt','gene','description','effect','impact','z.fsts.lc.se','z.fsts.lc.lh','z.fsts.se.lh')]
colnames(s2) <- c('chrom','pos','ref','alt','gene','description','effect','impact','z.fst.lc.se','z.fst.lc.lh','z.fst.se.lh')
s2 <- s2[order(s2$chrom, as.numeric(s2$pos)),]
write.table(s2, '/Users/Avril/Desktop/S2.txt', row.names=FALSE, quote=FALSE, sep='\t')

## S3: SnpEff annotations for hatchery-associated outlier SNPs
## Columns: chrom, pos, ref, alt, gene, description, effect, impact, z.fst.lc.se, z.fst.lc.lh, z.fst.se.lh
s3 <- trans.snps[which(trans.snps$loc %in% hilite),]
s3 <- s3[which(s3$z.fsts.lc.se >= 3 & s3$z.fsts.lc.lh >= 3),]
s3 <- s3[,c('chrom','pos','ref','alt','gene','description','effect','impact','z.fsts.lc.se','z.fsts.lc.lh','z.fsts.se.lh')]
colnames(s3) <- c('chrom','pos','ref','alt','gene','description','effect','impact','z.fst.lc.se','z.fst.lc.lh','z.fst.se.lh')
s3 <- s3[order(s3$chrom, as.numeric(s3$pos)),]
write.table(s3, '/Users/Avril/Desktop/S3.txt', row.names=FALSE, quote=FALSE, sep='\t')

## S4: SnpEff annotations for life history-associated outlier SNPs
## Columns: chrom, pos, ref, alt, gene, description, effect, impact, z.fst.lc.se, z.fst.lc.lh, z.fst.se.lh
s4 <- trans.snps[which(trans.snps$loc %in% hilite),]
s4 <- s4[which(s4$z.fsts.se.lh >= 3 & s4$z.fsts.lc.lh >= 3),]
s4 <- s4[,c('chrom','pos','ref','alt','gene','description','effect','impact','z.fsts.lc.se','z.fsts.lc.lh','z.fsts.se.lh')]
colnames(s4) <- c('chrom','pos','ref','alt','gene','description','effect','impact','z.fst.lc.se','z.fst.lc.lh','z.fst.se.lh')
s4 <- s4[order(s4$chrom, as.numeric(s4$pos)),]
write.table(s4, '/Users/Avril/Desktop/S4.txt', row.names=FALSE, quote=FALSE, sep='\t')
