#PBS -N download_sra
#PBS -q beagle
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=300:00:00

module purge
module load sra-toolkit/2.9.2

cd /scratch/snyder/h/harder/pop_comp/sra_downloads/ #$PBS_O_WORKDIR

cat ../scripts/sra_access_list.txt | while read line
	do
	fastq-dump --split-files --origfmt --gzip ${line}
done
