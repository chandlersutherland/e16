#!/bin/bash
#SBATCH --job-name=fastq_download.mo17
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

## download the sra files for epigenomes, rnaseq, and atacseq 
## B73 
sample='MO17'

module load python
module load sra-tools

SCRATCH_DIR='/global/scratch/users/chandlersutherland/e16'

atac='SRR13920262 SRR13920263'
bisulfite='SRR12463980 SRR12463981'
rna='SRR12454910 SRR12454911 SRR12454912 SRR12454913 SRR12454915'

atac_output=${SCRATCH_DIR}/${sample}/atac_fastq_files
mkdir -p $atac_output

#for i in $atac
#	do 
#	fasterq-dump --threads $SLURM_NTASKS -O $atac_output -t $SCRATCH_DIR -p $i
#	echo "$i bisulfite for sample $sample complete"
#done 

bs_output=${SCRATCH_DIR}/${sample}/bs_fastq_files
mkdir -p $bs_output

for i in $bisulfite
	do 
	fasterq-dump --threads $SLURM_NTASKS -O $bs_output -t $SCRATCH_DIR -p $i
	echo "$i bisulfite for sample $sample complete"
done 

rna_output=${SCRATCH_DIR}/${sample}/rna_fastq_files
mkdir -p $rna_output

for i in $rna
	do 
	fasterq-dump --threads $SLURM_NTASKS -O $rna_output -t $SCRATCH_DIR -p $i
	echo "$i rna for sample $sample complete"
done 

