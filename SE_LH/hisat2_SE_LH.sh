#!/bin/bash

#PBS -N hisat2+samtools
#PBS -q beagle
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

cd $PBS_O_WORKDIR

module purge
module load hisat2/2.1.0
module load samtools/1.8


cd /scratch/snyder/h/harder/pop_comp/hisat2/ #$PBS_O_WORKDIR

cat ../scripts/sra_access_list.txt | while read line
	do

	hisat2 --dta --time --threads 20 \
	-S ${line}.sam \
	-x simple_names_chroms_scaf_MT \
	-1 ../sra_downloads/${line}_1.fastq.gz \
	-2 ../sra_downloads/${line}_2.fastq.gz

	samtools view -@ 20 -Su ${line}.sam | samtools sort -@ 20 -o ${line}.sorted.bam

done
