#!/bin/bash

#PBS -N genomics_dbi
#PBS -q darwin
#PBS -l nodes=1:ppn=10,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

cd $PBS_O_WORKDIR

module purge
module load GATK/4

ref=/scratch/snyder/h/harder/ssalar_genome_and_annotation_files/genome_and_trans_and_indices/simple_names_chroms_scaf_MT.fa

cat chroms_only.txt | while read line 
	do 
	gatk GenomicsDBImport \
		-V split_deduped_wRG__1_hisat2_1b.sorted.bam.g.vcf \
		-V split_deduped_wRG__1_hisat2_6-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__1_hisat2_8-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__1_hisat2_11-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__2_hisat2_1c.sorted.bam.g.vcf \
		-V split_deduped_wRG__2_hisat2_6a.sorted.bam.g.vcf \
		-V split_deduped_wRG__2_hisat2_8b.sorted.bam.g.vcf \
		-V split_deduped_wRG__2_hisat2_11a.sorted.bam.g.vcf \
		-V split_deduped_wRG__3_hisat2_3-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__3_hisat2_6-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__3_hisat2_9-a2.sorted.bam.g.vcf \
		-V split_deduped_wRG__3_hisat2_11-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__4_hisat2_3a.sorted.bam.g.vcf \
		-V split_deduped_wRG__4_hisat2_6b.sorted.bam.g.vcf \
		-V split_deduped_wRG__4_hisat2_9a.sorted.bam.g.vcf \
		-V split_deduped_wRG__4_hisat2_11b.sorted.bam.g.vcf \
		-V split_deduped_wRG__5_hisat2_3-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__5_hisat2_7-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__5_hisat2_9-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__5_hisat2_13a2.sorted.bam.g.vcf \
		-V split_deduped_wRG__6_hisat2_3b.sorted.bam.g.vcf \
		-V split_deduped_wRG__6_hisat2_7a.sorted.bam.g.vcf \
		-V split_deduped_wRG__6_hisat2_9b.sorted.bam.g.vcf \
		-V split_deduped_wRG__6_hisat2_13-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__7_hisat2_4a2.sorted.bam.g.vcf \
		-V split_deduped_wRG__7_hisat2_7-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__7_hisat2_13-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__8_hisat2_4-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__8_hisat2_7b.sorted.bam.g.vcf \
		-V split_deduped_wRG__8_hisat2_13b.sorted.bam.g.vcf \
		-V split_deduped_wRG__9_hisat2_1-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__9_hisat2_4-b2.sorted.bam.g.vcf \
		-V split_deduped_wRG__9_hisat2_8-a.sorted.bam.g.vcf \
		-V split_deduped_wRG__10_hisat2_1-b.sorted.bam.g.vcf \
		-V split_deduped_wRG__10_hisat2_4b.sorted.bam.g.vcf \
		-V split_deduped_wRG__10_hisat2_8a.sorted.bam.g.vcf \
		--genomicsdb-workspace-path "$line"_database_gatk_genomics \
		--intervals $line
	
	gatk GenotypeGVCFs \
		-R $ref \
		-V gendb://"$line"_database_gatk_genomics \
		-new-qual true \
		-O "$line"_genotype_output.vcf
		
	done
