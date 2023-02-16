#!/bin/bash
#SBATCH --job-name=fastq_download.b73
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=09:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

## download the sra files for epigenomes, rnaseq, and atacseq 
## B73 
sample='B73'

module load python
module load sra-tools

SCRATCH_DIR='/global/scratch/users/chandlersutherland/e16'

atac='SRR13920264 SRR13920265'
bisulfite='SRR12463972 SRR12463973'
rna='SRR12454914 SRR12454925 SRR12454936 SRR12454937'

#atac_output=${SCRATCH_DIR}/${sample}/atac_fastq_files
#mkdir -p $atac_output

#for i in $atac
#	do 
#	fasterq-dump --threads $SLURM_NTASKS -O $atac_output -t $SCRATCH_DIR -p $i
#	echo "$i atac for sample $sample complete"
#done 

bs_output=${SCRATCH_DIR}/${sample}/bs_fastq_files
mkdir -p $bs_output

#for i in $bisulfite
#	do 
#	fasterq-dump --threads $SLURM_NTASKS -O $bs_output -t $SCRATCH_DIR -p $i
#	echo "$i bisulfite for sample $sample complete"
#done 

time fasterq-dump --threads $SLURM_NTASKS -O $bs_output -t $SCRATCH_DIR -p SRR12463973
echo "SRR12463973 for sample $sample complete"

rna_output=${SCRATCH_DIR}/${sample}/rna_fastq_files
mkdir -p $rna_output

for i in $rna
	do 
	time fasterq-dump --threads $SLURM_NTASKS -O $rna_output -t $SCRATCH_DIR -p $i
	echo "$i rna for sample $sample complete"
done 

