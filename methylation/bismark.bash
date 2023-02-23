#!/bin/bash
#SBATCH --job-name=bismark
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load bowtie2
module load Bismark
module load samtools
module load parallel 

trim_output=$base/$sample/em/trimmed
genome=$base/$sample/genome/bismark
bismark_output=$base/$sample/em/bismark
mkdir -p $bismark_output
#need to export three variables, trim_output directory, which is our input directory, bismark_output directory, and the genome file 
cd $trim_output

a=$(find . -type f -name '*_1_val_1.fq')
accession=$(basename -s _1_val_1.fq $a)

BISMARK() {
    echo "beginning bismark on sample ${sample} ${1}"
    bismark --genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	-p 4 \
	-1 "${1}"_1_val_1.fq -2 "${1}"_2_val_2.fq
}

export genome=$genome
export bismark_output=$bismark_output
export -f BISMARK

parallel BISMARK ::: $accession
