#!/bin/bash

#PBS -N gatkbp_HISAT2_group_
#PBS -q darwin
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

cd $PBS_O_WORKDIR

module purge

module load bioinfo
module load picard-tools/2.18.2
module load samtools/1.8
module load GATK/4

##----------------------------------------------------------------------------------------
## add read group information, as required by GATK for downstream analysis

cd /scratch/snyder/h/harder/analyses_ssalar_rna_seq/HISAT2_gatk_bp/sorted_hisat2_bams

for i in *_group_*.sorted.bam
 	do
 	gatk AddOrReplaceReadGroups -I $i -O /scratch/snyder/h/harder/analyses_ssalar_rna_seq/HISAT2_gatk_bp/_group_/wRG_$i -LB $i -PL illumina -PU A -SM $i
 	done

cd /scratch/snyder/h/harder/analyses_ssalar_rna_seq/HISAT2_gatk_bp/_group_/
 
##----------------------------------------------------------------------------------------
## mark duplicates 

for i in wRG*.bam
	do
	gatk MarkDuplicates -I $i -O deduped_$i -M $i.metrics
	done
	
##----------------------------------------------------------------------------------------
## create index

for i in deduped*.bam
	do
	gatk BuildBamIndex -I $i
	done

##----------------------------------------------------------------------------------------
## split & trim reads = splits reads into exons and hard-clips any sequences overhanging
## into the intronic regions

module unload GATK/4
module load GATK/3.8.0

for i in deduped*.bam
	do
	GenomeAnalysisTK -T SplitNCigarReads -R /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/genome_and_trans_and_indices/simple_names_chroms_scaf_MT.fa -nt 1 -I $i -o split_$i -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
	done

##----------------------------------------------------------------------------------------
## call variants, reporting all genotypes for later joint genotyping

module unload GATK/3.8.0
module load GATK/4

for i in split*.bam
	do
	gatk HaplotypeCaller -R /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/genome_and_trans_and_indices/simple_names_chroms_scaf_MT.fa -I $i --dont-use-soft-clipped-bases --emit-ref-confidence GVCF -stand-call-conf 20.0 -O $i.g.vcf
	done

##----------------------------------------------------------------------------------------
