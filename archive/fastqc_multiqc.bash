#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load fastqc 
source activate e16 

SCRATCH_DIR=/global/scratch/users/chandlersutherland/e16
OUTPUT_DIR=$SCRATCH_DIR/fastqc

mkdir -p $OUTPUT_DIR

cd $SCRATCH_DIR
FILES=$(find . -type f -name '*fastq' -print)

fastqc -o $OUTPUT_DIR -t 24 $FILES 

multiqc -f $OUTPUT_DIR
