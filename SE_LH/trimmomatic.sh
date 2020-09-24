#PBS -N trimmomatic+qc
#PBS -q beagle
#PBS -l nodes=1:ppn=20,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

cd $PBS_O_WORKDIR

module purge
module load bioinfo
module load trimmomatic/0.36
module load fastqc

cd /scratch/snyder/h/harder/pop_comp/sra_downloads/ #$PBS_O_WORKDIR

mkdir fastqc_before
fastqc -t 20 -o fastqc_before *fastq.gz

trimmer="LEADING:20 TRAILING:20 MINLEN:30"

cat ../scripts/sra_access_list.txt | while read line
	do
	trimmomatic PE -threads 20 \
	${line}_1.fastq.gz ${line}_2.fastq.gz \
	${line}_1.paired.filtered.fastq.gz ${line}_1.unpaired.filtered.fastq.gz \
	${line}_2.paired.filtered.fastq.gz ${line}_2.unpaired.filtered.fastq.gz \
	$trimmer
done

mkdir fastqc_after
fastqc -t 20 -a ../scripts/adapters.fa -o fastqc_after *.filtered.fastq.gz
