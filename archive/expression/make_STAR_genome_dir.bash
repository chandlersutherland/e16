#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load pigz 
module load gffread 

#first, unzip all the genomes. This doens't seem to have happened during bismark genome prep... 
GENOME=$base/$sample/genome
cd $GENOME
#unpigz *.gz 

gff3=$(find . -type f -name '*.gff3' -print)
echo $gff3
for f in $gff3 
do 
	basename=$(basename -s .gff3 $f)
	gffread -T -v -o $basename.gtf $f
done 

echo 'converted $sample annotations to gtf format' 

source activate e16
#the genome directory is the output of the STAR converted genome index 
GENOME_DIR=$base/$sample/genome/STAR
mkdir -p $GENOME_DIR

#name the input file variables 
FASTA_FILE=$GENOME/*.fa
GTF_FILE=$GENOME/*.1.gtf

#The sjdb overhang depends on the read length. Here we have 150bp reads, so set to 149 
STAR --runThreadN $SLURM_NTASKS \
	 --runMode genomeGenerate \
	 --genomeDir $GENOME_DIR \
	 --genomeFastaFiles $FASTA_FILE \
	 --sjdbGTFfile $GTF_FILE  \
	 --sjdbOverhang 149

echo 'STAR genome finished for $sample' 
