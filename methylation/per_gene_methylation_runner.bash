#!/bin/bash
#SBATCH --job-name=per_gene_meth_runner
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#this script passes a bam file to the python script samtools_coverage, which outputs a file $BASENAME_clean_coverage.tsv which gives the mean depth over the NLRs
#can easily be made into a for loop to pass multiple files through 
#pass variable coverage_input, which should have a directory of sam files ready for coverage processing 

module load python

#point to the all gene bed file 
gene_positions=/global/scratch/users/chandlersutherland/e16/${sample}/genome/*_all_gene.bed

#where are the coverage files that have not been processed yet
cov_dir=/global/scratch/users/chandlersutherland/e16/${sample}/em/bedGraph_highcov

cd $cov_dir
#First define the finished files
finished=$(find . -type f -name '*.tsv')
prefix=$(basename -s "_per_gene_met_${sample}.tsv" $finished)
for i in $prefix; do echo $i >> finished_prefix; done 



#get all of the files 
all=$(find . -type f -name '*.cov')
all_prefix=$(basename -s '.bed.gz.bismark.cov' $all)
for i in $all_prefix; do echo $i >> all_prefix; done 

#define the complement, aka all unknown files 
unfinished=$(comm -23 <(sort all_prefix) <(sort finished_prefix))

echo "$unfinished have not been completed. Starting meth extraction now"
rm prefix 
rm finished_prefix 

for f in $unfinished
do 
	python $HOME/e16/methylation/per_gene_methylation.py $cov_dir/$f.bed.gz.bismark.cov $gene_positions 
done
