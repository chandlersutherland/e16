#!/bin/bash
#SBATCH --job-name=expression_processing
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

cd /global/home/users/chandlersutherland/e16/
base="/global/scratch/users/chandlersutherland/e16"
module load python 

#first, generate a STAR genome 
#while read sample; do sbatch --job-name=$sample.genome --export=base=$base,sample=$sample \
#-A fc_kvkallow expression/make_STAR_genome_dir.bash; done < sample.txt 
sample='CML103' 

#unpigz files 
while read tissue; do sbatch --job-name=$sample.$tissue.unpigz \
--export=base=$base,sample=$sample,tissue=$tissue -A co_minium -p savio4_htc \
--qos minium_htc4_normal unpigz.sh; done < tissues.txt

sleep 5m

#run STAR in quant mode 
while read tissue; do sbatch --job-name=$sample.$tissue.STAR --export=base=$base,sample=$sample,tissue=$tissue \
-A co_minium -p savio4_htc --qos minium_htc4_normal expression/STAR.bash; done < tissues.txt 

sleep 15m
#create output file 
python tpm_calc.py $sample 

#QC coverage of each NLR 
#while read sample; do sbatch --job-name=$sample.coverage --export=base=$base,sample=$sample \
#-A co_minium --partition savio4_htc expression/samtools_coverage_runner.bash; done < sample.txt 