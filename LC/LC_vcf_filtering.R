## script for further filtering LC SNPs and calling genotypes

setwd('/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/snps_and_gene_expression/variant_calling/20_12_15_testing_scripts_for_ms_upload/')

library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)

`%notin%` <- Negate(`%in%`)

##### LOAD AND INSPECT DATA #####
### read in variant file produced by variant_filtration.sh 
vcf <- read.vcfR('/Users/Avril/Desktop/.scratch/analyses_ssalar_rna_seq/HISAT2_gatk_bp/stacks/gvcf_files/all_chroms_all_samps.vcf.gz')

vcf@fix[1:10,1:ncol(vcf@fix)]
vcf@fix[1:10,1:ncol(vcf@fix)-1]
nrow(vcf@fix) ## 2,494,061 SNPs

## check filter results for all SNPs
## FS = Fisher Strand values > 30.0; higher values indicate strand bias more likely, indicative of false positive calls
## QD = Qual By Depth values < 2.0; variant call confidence normalized by depth of reads supporting a variant
## SnpCluster = cluster of at least 3 SNPs within a window of 35 bases between them
## PASS = no filters violated
table(vcf@fix[,7]) ## 1,948,479 PASS
vcf
vcf <- vcf[which(vcf@fix[,7]=="PASS"),] ## retain only PASS SNPs
vcf

## extract indels - removes any line w/ an allele (REF or ALT) > 1 character in length
vcf <- extract.indels(vcf, return.indels=FALSE)
vcf ## 1,652,102 

## keep biallelic SNPs only
vcf <- vcf[is.biallelic(vcf),] 
vcf ## 1,642,500 

## check allele depths per sample w/ plot
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)

##### STRINGENT DEPTH FILTERING APPROACH USING MIN COUNTS, rowSums(), and rowMeans() #####
n <- 18

## keep SNPs where n samples have DP > 5
dim(dp[rowSums(dp > 5) >= n,]) 
dp <- dp[rowSums(dp > 5) >= n,] 
dim(dp) 

vcf@fix <- cbind(vcf@fix, paste0(vcf@fix[,1],'_',vcf@fix[,2]))
vcf <- vcf[which(vcf@fix[,9] %in% rownames(dp)),]
vcf@fix <- vcf@fix[,-c(9,10)]

vcf
gtypes <- extract.gt(vcf, element='GT') ## table(gtypes[,c(1:36)]) --> all samples gtyped at all loci
gtypes <- gsub("0/0", "0", gtypes, fixed=TRUE)
gtypes <- gsub("0|0", "0", gtypes, fixed=TRUE)
gtypes <- gsub("0/1", "1", gtypes, fixed=TRUE)
gtypes <- gsub("0|1", "1", gtypes, fixed=TRUE)
gtypes <- gsub("1/0", "1", gtypes, fixed=TRUE)
gtypes <- gsub("1|0", "1", gtypes, fixed=TRUE)
gtypes <- gsub("1/1", "2", gtypes, fixed=TRUE)
gtypes <- gsub("1|1", "2", gtypes, fixed=TRUE)
colnames(gtypes) <- c('1-b','4b','8a','11-a','1b','6-a','8-b','11a','1c','6a','8b','11-b','3-a','6-b','9-a2','11b','3a','6b','9a','13a2','3-b','7-a','9-b','13-a','3b','7a','9b','13-b','4a2','7-b','13b','4-a','7b','1-a','4-b2','8-a')

names <- do.call(rbind, strsplit(rownames(gtypes), split='_', fixed=TRUE))
rownames(gtypes) <- paste0(names[,1],'_',names[,2],'|',names[,3])

##### !! do some gtype filtering based on individual read depth #####
## Brouard et al. 2020: require depth=5 for genotype calls
colnames(dp) <- c('1-b','4b','8a','11-a','1b','6-a','8-b','11a','1c','6a','8b','11-b','3-a','6-b','9-a2','11b','3a','6b','9a','13a2','3-b','7-a','9-b','13-a','3b','7a','9b','13-b','4a2','7-b','13b','4-a','7b','1-a','4-b2','8-a')
## limit down to genotyped loci
rownames(dp) <- gsub('.1_','.1|', rownames(dp))
temp.dp <- dp[which(rownames(dp) %in% rownames(gtypes)),]
all(rownames(temp.dp) == rownames(gtypes)) ## TRUE
all(colnames(temp.dp) == colnames(gtypes)) ## TRUE

## identify gtype calls supported with <5 reads and set gtypes to NA
temp.gtypes <- gtypes
for(i in 1:ncol(temp.dp)){
  dps <- temp.dp[,i]
  loci <- names(dps[which(dps < 5)])
  temp.gtypes[which(rownames(temp.gtypes) %in% loci), i] <- 'NA'
  print(i)
}
## limit loci to those with gtypes in at least 20 samples
num.nas <- rowSums(is.na(temp.gtypes))
num.nas <- num.nas[which(num.nas >= 16)]
temp.gtypes <- temp.gtypes[which(rownames(temp.gtypes) %notin% names(num.nas)),]

## write new gtype file with samples/loci with RD < 5 set to NA gtype
write.csv(temp.gtypes, '20_12_28_genotype_data_36_samples.csv', quote=FALSE) ## 234,012 SNPs









