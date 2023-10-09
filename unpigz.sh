#!/bin/bash
#SBATCH --job-name=unpigz
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:10:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=FAIL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load pigz
INPUT=$base/$sample/rna_${tissue}

cd $INPUT

files=$(find . -type f -name '*.gz')
unpigz $files 
