setwd('/Users/Avril/Desktop/all_ch_3_dirs/02_output_tables_all_analyses/')

library(rtracklayer)
library(vcfR)
library(scales)
library(ggrepel)
library(dplyr)
library(ade4)
library(adegenet)
library(VennDiagram)
library(TeachingDemos)

load('combining_LC_and_SE_LH_snps_thru_kNN.RData') ## or run 3_population_analyses.R first

## set gene of interest
gene <- 'LOC106610379' ## zinc finger and BTB domain-containing protein 18-like, gene with 4 kNN missense outlier SNPs
# gene <- 'LOC106576177' ## socs2
# gene <- 'LOC106571946' ## AC6
## get reference info
temp.ass <- sub.ass[which(sub.ass$gene == gene),]
ref.gene <- temp.ass[which(temp.ass$type == 'gene'),]
ref.exons <- temp.ass[which(temp.ass$type == 'exon'),]
## get SNPs in and near gene (5kb up- and downstream)
temp.snps <- pop.afs[which(pop.afs$chrom %in% temp.ass$seqid),]
temp.snps$pos <- as.numeric(temp.snps$pos)
temp.snps <- temp.snps[which(temp.snps$pos < (ref.gene$end + 5000)),]
temp.snps <- temp.snps[which(temp.snps$pos > (ref.gene$start - 5000)),]
knn.out <- knn.out[which(knn.out$zfst.out == TRUE),]
temp.out <- temp.snps[which(temp.snps$loc %in% knn.out$loc),]
temp.not.out <- temp.snps[which(temp.snps$loc %notin% temp.out$loc),]
msense.snps <- trans.snps[which(trans.snps$loc %in% temp.out$loc & trans.snps$effect == 'missense_variant'), 'loc']

## plot allele frequencies along the gene
pdf('/Users/Avril/Desktop/gene_diagram.pdf', width=7, height=4)
par(mar=c(5.1, 5.1, 4.1, 2.1))
# plot(temp.snps$pos, temp.snps$lc.ref.af, xlim=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100), ylim=c(0,1.36), ## limited by SNP position -- or --
plot(temp.snps$pos, temp.snps$lc.ref.af, xlim=c(ref.gene$start, ref.gene$end), ylim=c(0,1.36), ## show whole gene
     pch=19, col='white', ylab='Reference allele\nfrequency', xlab='Chromosome position (Mbp)', main=gene, yaxt='n', bty='n', xaxt='n')
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(-0.5, -0.5, 1.3, 1.3), col='gray96', border=NA)
  }
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(1.3, 1.3, 1.35, 1.35), col='black', border=NA)
  }
  mtext(text='exons', side=2, at=1.325, line=1, las=2)
  axis(1, at=c(ref.gene$start, ref.gene$end), ## show whole gene
       labels=c(round(ref.gene$start/1e6, digits=4),
                round(ref.gene$end/1e6, digits=4)))
  axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c('0.0', NA, '0.5', NA, '1.0'))
  ## plot non-outlier SNPs
  points(temp.not.out$pos, temp.not.out$lc.ref.af, col=alpha(lc.col, 0.8), pch=21, bg='white', cex=1.1, lwd=3)
  points(temp.not.out$pos, temp.not.out$se.ref.af, col=alpha(se.col, 0.8), pch=21, bg='white', cex=1.1, lwd=3)
  points(temp.not.out$pos, temp.not.out$lh.ref.af, col=alpha(lh.col, 0.8), pch=21, bg='white', cex=1.1, lwd=3)
  ## plot outlier SN
  points(temp.out$pos, temp.out$lc.ref.af, col=alpha(lc.col, 0.8), pch=19, cex=1.1)
  points(temp.out$pos, temp.out$se.ref.af, col=alpha(se.col, 0.8), pch=19, cex=1.1)
  points(temp.out$pos, temp.out$lh.ref.af, col=alpha(lh.col, 0.8), pch=19, cex=1.1)
  

## plot legend separately
plot(temp.snps$pos, temp.snps$lc.ref.af, xlim=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100), ylim=c(0,1.36),
       pch=19, col='white', ylab='Reference allele\nfrequency', xlab='Chromosome position (Mbp)', main=gene, yaxt='n', bty='n', xaxt='n')
  # legend(x=mean(min(temp.snps$pos), max(temp.snps$pos)), y=1, 
  #        legend=c('LC','SE','LH','Outlier SNP'), pch=c(19,19,19,NA), col=c(lc.col, se.col, lh.col, 'grey30'),
  #        lty=c(NA,NA,NA,2), lwd=c(NA,NA,NA,0.5), pt.cex=c(1.1,1.1,1.1,NA), box.lty=1, box.lwd=0.5, bg='transparent')
  legend(x=mean(min(temp.snps$pos), max(temp.snps$pos)), y=1, 
         legend=c('LC','SE','LH'), pch=c(19,19,19), col=c(lc.col, se.col, lh.col), pt.cex=c(1.1,1.1,1.1), box.lty=1, box.lwd=0.5, bg='transparent')
dev.off()
    
## plot AF for each SNP
pdf('/Users/Avril/Desktop/pop_af.pdf', width=4, height=3)
par(mar=c(4.1, 4.5, 1.1, 1.1))
plot(x=c(rep(1, nrow(temp.snps)), rep(2, nrow(temp.snps)), rep(3, nrow(temp.snps))),
     y=c(temp.snps$lc.ref.af, temp.snps$se.ref.af, temp.snps$lh.ref.af), pch=19, ylim=c(0,1), yaxt='n',
     col='white', xlab='Population', ylab='Reference allele frequency', xaxt='n')
  axis(1, at=c(1,2,3), labels=c('LC', 'SE', 'LH'))
  axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c('0.0', NA, '0.5', NA, '1.0'))
  for(i in 1:nrow(temp.snps)){
    if(temp.snps$loc[i] %in% msense.snps){
      lines(x=c(1,2,3), y=c(temp.snps$lc.ref.af[i], temp.snps$se.ref.af[i], temp.snps$lh.ref.af[i]),
          col='red', lwd=1)
    }
  }
  for(i in 1:nrow(temp.snps)){
    if(temp.snps$loc[i] %notin% msense.snps){
      lines(x=c(1,2,3), y=c(temp.snps$lc.ref.af[i], temp.snps$se.ref.af[i], temp.snps$lh.ref.af[i]),
          col='black', lty=2, lwd=1)
    }
  }    
  points(x=c(rep(1, nrow(temp.snps)), rep(2, nrow(temp.snps)), rep(3, nrow(temp.snps))),
         y=c(temp.snps$lc.ref.af, temp.snps$se.ref.af, temp.snps$lh.ref.af), pch=19, cex=1.5,
         col=c(rep(lc.col, nrow(temp.snps)), rep(se.col, nrow(temp.snps)), rep(lh.col, nrow(temp.snps))))
dev.off()

## plot Fst for each SNP
pdf('/Users/Avril/Desktop/gene_Fsts.pdf', width=4, height=3.5)
par(mar=c(4.1, 4.5, 0, 1.1))
## LC vs LH
# plot(temp.snps$pos, temp.snps$fst.lc.lh, xlim=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100), ylim=c(-0.05,1.36), ## limited by SNP position -- or --
plot(temp.snps$pos, temp.snps$fst.lc.lh, xlim=c(ref.gene$start, ref.gene$end), ylim=c(-0.05,1.36), ## show whole gene
     pch=19, col='white', ylab=bquote(italic(F)[ST]), xlab='Chromosome position (Mbp)', yaxt='n', bty='n', xaxt='n')
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(-0.5, -0.5, 1.3, 1.3), col='gray96', border=NA)
  }
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(1.3, 1.3, 1.35, 1.35), col='black', border='black')
  }
  cds <- temp.ass[which(temp.ass$type == 'CDS'),]
  for(i in 1:nrow(cds)){
    polygon(x=c(cds$start[i], cds$end[i], cds$end[i], cds$start[i]),
            y=c(1.3, 1.3, 1.35, 1.35), col='white', border='black')
  }
  mtext(text='exons', side=2, at=1.325, line=1, las=2)
  # axis(1, at=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100), ## limited by SNP position --or-- 
  #      labels=c(round((min(temp.snps$pos) - 100)/1e6, digits=4),
  #               round((max(temp.snps$pos) + 100)/1e6, digits=4)))
  axis(1, at=c(ref.gene$start, ref.gene$end), ## show whole gene
       labels=c(round(ref.gene$start/1e6, digits=4),
                round(ref.gene$end/1e6, digits=4)))
  axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c('0.0', NA, '0.5', NA, '1.0'))
  ## plot non-outlier SNPs
  points(temp.out$pos, temp.out$fst.lc.lh, col=alpha(lc.vs.lh, 0.8), pch=21, bg='white', cex=1.1, lwd=3)
  ## plot outlier SN
  points(temp.not.out$pos, temp.not.out$fst.lc.lh, col=alpha(lc.vs.lh, 0.8), pch=19, cex=1.1)
  
  ## lines for missense SNPs
  for(j in msense.snps){
    sub <- temp.out[which(temp.out$loc ==j),]
    lines(x=c(sub$pos, sub$pos), y=c(1.15, 1.25), lwd=3)
  }
  mtext(text='missense', side=2, at=1.2, line=1, las=2)
## plot legend separately
plot(temp.snps$pos, temp.snps$fst.lc.lh, xlim=c(ref.gene$start, ref.gene$end), ylim=c(-0.05,1.36), ## show whole gene
     pch=19, col='transparent', ylab=bquote(italic(F[ST])[(LC-LH)]), xlab='Chromosome position (Mbp)', yaxt='n', bty='n', xaxt='n')
  legend('center', legend=c('CDS','UTR','outlier'), fill=c('white','black',NA), pch=c(NA,NA,21), pt.bg=c(NA,NA,'white'), col=c(NA,NA,'black'), pt.lwd=c(NA,NA,3), pt.cex=c(NA,NA,1.1))
## plot legend separately
plot(temp.snps$pos, temp.snps$fst.lc.lh, xlim=c(ref.gene$start, ref.gene$end), ylim=c(-0.05,1.36), ## show whole gene
       pch=19, col='transparent', ylab=bquote(italic(F[ST])[(LC-LH)]), xlab='Chromosome position (Mbp)', yaxt='n', bty='n', xaxt='n')
  legend('center', legend=c('outlier'), pch=c(21), pt.bg=c('white'), col=c('black'), pt.lwd=c(NA,NA,3), pt.cex=c(NA,NA,1.1))
  
## LC vs SE
plot(temp.snps$pos, temp.snps$fst.lc.se, xlim=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100), ylim=c(-0.05,1.36),
       pch=19, col='white', ylab=bquote(italic(F)[ST]), xlab='Chromosome position (Mbp)', main=gene, yaxt='n', bty='n', xaxt='n')
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(-0.5, -0.5, 1.2, 1.2), col='gray96', border=NA)
  }
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(1.2, 1.2, 1.25, 1.25), col='black', border=NA)
  }
  axis(1, at=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100),
       labels=c(round((min(temp.snps$pos) - 100)/1e6, digits=4),
                round((max(temp.snps$pos) + 100)/1e6, digits=4)))
  axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c('0.0', NA, '0.5', NA, '1.0'))
  ## plot non-outlier SNPs
  points(temp.out$pos, temp.out$fst.lc.se, col=alpha(lc.vs.se, 0.8), pch=21, bg='white', cex=1.1, lwd=3)
  ## plot outlier SN
  points(temp.not.out$pos, temp.not.out$fst.lc.se, col=alpha(lc.vs.se, 0.8), pch=19, cex=1.1)
  mtext(text='exons', side=2, at=1.225, line=1, las=2)
## SE vs LH
plot(temp.snps$pos, temp.snps$fst.lc.lh, xlim=c(ref.gene$start, ref.gene$end), ylim=c(-0.05,1.36), ## show whole gene
# plot(temp.snps$pos, temp.snps$fst.se.lh, xlim=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100), ylim=c(-0.05,1.36),
       pch=19, col='white', ylab=bquote(italic(F[ST])), xlab='Chromosome position (Mbp)', main=gene, yaxt='n', bty='n', xaxt='n')
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(-0.5, -0.5, 1.2, 1.2), col='gray96', border=NA)
  }
  for(i in 1:nrow(ref.exons)){
    polygon(x=c(ref.exons$start[i], ref.exons$end[i], ref.exons$end[i], ref.exons$start[i]),
            y=c(1.2, 1.2, 1.25, 1.25), col='black', border=NA)
  }
  # axis(1, at=c(min(temp.snps$pos) - 100, max(temp.snps$pos) + 100), ## limited by SNP position --or-- 
  #      labels=c(round((min(temp.snps$pos) - 100)/1e6, digits=4),
  #               round((max(temp.snps$pos) + 100)/1e6, digits=4)))
  axis(1, at=c(ref.gene$start, ref.gene$end), ## show whole gene
       labels=c(round(ref.gene$start/1e6, digits=4),
                round(ref.gene$end/1e6, digits=4)))
  axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c('0.0', NA, '0.5', NA, '1.0'))
  ## plot non-outlier SNPs
  points(temp.out$pos, temp.out$fst.se.lh, col=alpha(se.vs.lh, 0.8), pch=21, bg='white', cex=1.1, lwd=3)
  ## plot outlier SN
  points(temp.not.out$pos, temp.not.out$fst.se.lh, col=alpha(se.vs.lh, 0.8), pch=19, cex=1.1)
  mtext(text='exons', side=2, at=1.225, line=1, las=2)
dev.off()

  
