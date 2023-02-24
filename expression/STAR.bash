#!/bin/bash
#SBATCH --job-name=star
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=04:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python
module load parallel
source activate e16

#define input variables: input and output directories, number of threads 
INPUT=$base/$sample/rna_tip
cd $INPUT
a=$(find . -type f -name '*_1.fastq')
accession=$(basename -s _1.fastq $a)
reps=$(echo "$accession" | wc -l)
threads=$(echo `expr $SLURM_NTASKS / $reps`)

GENOME_DIR=$base/$sample/$genome/$STAR
STAR_OUTPUT=$INPUT/STAR
mkdir -p $STAR_OUTPUT

#define star function, which takes in a accession name and writes the STAR output to the STAR_OUTPUT directory 
STAR_RUN (){ 
    STAR --runThreadN $threads \
		--genomeDir $GENOME_DIR \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts \
		--outFileNamePrefix "${STAR_OUTPUT}"/"${1}"_ \
		--readFilesIn "${1}"_1.fastq "${1}"_2.fastq 
    echo "finished $1"
}

#apply STAR_RUN to the input files
export GENOME_DIR=$GENOME_DIR
export STAR_OUTPUT=$STAR_OUTPUT
export threads=$threads
export -f STAR_RUN

parallel STAR_RUN ::: $accession