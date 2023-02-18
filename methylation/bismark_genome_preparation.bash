#!/bin/bash
#SBATCH --job-name=bismark_genome_preparation
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load bowtie2
module load samtools
module load Bismark 
module load python 
source activate e14 

#input is a genome folder, uncompress genome files 
GENOME=$base/$sample/genome
cd $GENOME
unpigz * 

#run genome preparation
time bismark_genome_preparation --verbose --bowtie2 --parallel 12 --path_to_aligner /global/home/groups/consultsw/sl-7.x86_64/modules/bowtie2/2.3.4.1 $GENOME
