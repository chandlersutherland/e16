#!/bin/bash
#SBATCH --job-name=fastq_download.w22
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
sample='W22'

module load python
module load sra-tools

SCRATCH_DIR='/global/scratch/users/chandlersutherland/e16'

atac='SRR13920266 SRR13920267'
bisulfite='SRR12463976 SRR12463977'
rna='SRR12454927 SRR12454928 SRR12454929 SRR12454930 SRR12454931'

atac_output=${SCRATCH_DIR}/${sample}/atac_fastq_files
mkdir -p $atac_output

#for i in $atac
#	do 
#	fasterq-dump --threads $SLURM_NTASKS -O $atac_output -t $SCRATCH_DIR -p $i
#	echo "$i atac for sample $sample complete"
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

