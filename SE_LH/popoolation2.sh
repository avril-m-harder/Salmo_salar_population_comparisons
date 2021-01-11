#!/bin/bash

#SBATCH --job-name=pop2_chroms
#SBATCH -A darwin
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 300:00:00

cd $SLURM_SUBMIT_DIR

module purge
module load bioinfo
module load popoolation2/1201
module load samtools/1.8

cd /scratch/snyder/h/harder/pop_comp/popoolation2_10_pools/

## index sorted BAM files
samtools index -@ 5 ./sorted_bams/SRR3986816.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986817.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986818.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986819.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986820.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986813.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986814.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986815.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986821.sorted.bam
samtools index -@ 5 ./sorted_bams/SRR3986822.sorted.bam

# ## run steps on individual chromosomes, leaving scaffolds out
cat chroms_only.txt | while read line

	do

 	## create mpileup file for each chromosome
 	samtools mpileup --ignore-RG -r ${line} \
 	--reference /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/genome_and_trans_and_indices/simple_names_chroms_scaf_MT.fa \
	./sorted_bams/SRR3986816.sorted.bam \
	./sorted_bams/SRR3986817.sorted.bam \
	./sorted_bams/SRR3986818.sorted.bam \
	./sorted_bams/SRR3986819.sorted.bam \
	./sorted_bams/SRR3986820.sorted.bam \
	./sorted_bams/SRR3986813.sorted.bam \
	./sorted_bams/SRR3986814.sorted.bam \
	./sorted_bams/SRR3986815.sorted.bam \
	./sorted_bams/SRR3986821.sorted.bam \
	./sorted_bams/SRR3986822.sorted.bam \
 	-o ./chrom_mpileups/"$line"_10_pools.mpileup

	##------ run PoPoolation2 -------
	## create synchronized file
	mpileup2sync.pl --input ./chrom_mpileups/"$line"_10_pools.mpileup --min-qual 20 \
	--fastq-type sanger --output ./chrom_syncs/"$line"_10_pools.sync

	done
	
cat ./chrom_syncs/*.sync >> all_chrom.sync

grep -Pv "0:0:0:0:0:0" all_chrom.sync > no_zeroes_all_chrom.sync