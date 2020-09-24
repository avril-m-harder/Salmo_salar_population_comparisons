#!/bin/bash

#PBS -N gatk_var_filt
#PBS -q beagle
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

cd $PBS_O_WORKDIR

module purge
module load GATK/4
module load bcftools/1.8

for i in *genotype_output.vcf
	do
	gatk VariantFiltration -R /scratch/snyder/h/harder/ssalar_genome_and_annotation_files/genome_and_trans_and_indices/simple_names_chroms_scaf_MT.fa -V $i -window 35 -cluster 3 --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O filtered_$i
	done

bcftools concat -f filtered_vcfs.txt -O z -o all_chroms_all_samps.vcf.gz --threads 20
