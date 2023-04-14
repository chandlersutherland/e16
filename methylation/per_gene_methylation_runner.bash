#!/bin/bash
#SBATCH --job-name=per_gene_meth_runner
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#this script passes a bam file to the python script samtools_coverage, which outputs a file $BASENAME_clean_coverage.tsv which gives the mean depth over the NLRs
#can easily be made into a for loop to pass multiple files through 
#pass variable coverage_input, which should have a directory of sam files ready for coverage processing 

module load python

#point to the all gene bed file 
gene_positions=/global/scratch/users/chandlersutherland/e16/${sample}/genome/*_all_gene.bed

#where are the coverage files 
cov_dir=/global/scratch/users/chandlersutherland/e16/${sample}/em/bedGraph_highcov
for f in $cov_dir/*.bed.gz.bismark.cov
do 
	python $HOME/nlr_features/methylation/per_gene_methylation.py $f $gene_positions 
done
