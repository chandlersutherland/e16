#!/bin/bash
#SBATCH --job-name=simulate_experiment_runner
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=03:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load r
module load r-packages

cd $SCRATCH

fasta=/global/scratch/users/chandlersutherland/e16/${sample}/genome/*.canonical.cds.fa
out_dir=/global/scratch/users/chandlersutherland/e16/${sample}/rna_tipe/simulated

mkdir -p $out_dir 

Rscript --no-save /global/home/users/chandlersutherland/e14/simulate_experiment.r $fasta $out_dir

#bash $HOME/e14/polyester_pipeline.bash 
