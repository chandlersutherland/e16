#!/bin/bash
#SBATCH --job-name=bismark_test
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=09:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load bowtie2
module load Bismark
module load samtools
module load parallel 

trim_output=$base/$sample/bs_fastq_files
genome=$base/$sample/genome
#need to export sample and base 
#FILES=$(find . -type f -name '*_1.fq' -print)
#ACCESSIONS=$(basename -a -s _1_val_1.fq $FILES)
#echo "beginning bismark on sample $sample"

bismark_output=$base/$sample/bismark_test

#just test one 
ACC='SRR12463972'

time bismark 
	--genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	--prefix nondirectional.trim1 \
	-p 4 \
	--non_directional \
	-1 $trim_output/trimmed/$ACC_1_val_1.fq -2 $trim_output/trimmed/$ACC_2_val_2.fq
	
time bismark 
	--genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	-p 4 \
	--prefix directional.trim1 \
	-1 $trim_output/trimmed/$ACC_1_val_1.fq -2 $trim_output/trimmed/$ACC_2_val_2.fq

time bismark 
	--genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	--prefix nondirectional.trim2 \
	-p 4 \
	--non_directional \
	-1 $trim_output/trimmed2/$ACC_1_val_1.fq -2 $trim_output/trimmed2/$ACC_2_val_2.fq
	
time bismark 
	--genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	-p 4 \
	--prefix directional.trim2 \
	-1 $trim_output/trimmed2/$ACC_1_val_1.fq -2 $trim_output/trimmed2/$ACC_2_val_2.fq
	
#export genome=$genome
#export bismark_output=$bismark_output
#export -f da_biz 

#make a final report 
cd $base/$sample/bismark_test
reports=$(find . -type f -name '*_PE_report.txt' -print)
for f in $reports 
	do 
		echo $f >> conglom_report.txt
		cat $f | grep -A 5 'Number of sequence pairs with unique' >> conglom_report.txt 
		cat $f | grep -A 5 'C methylated in CpG context:' >> conglom_report.txt
done 