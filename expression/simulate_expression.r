#!/usr/bin/env Rscript
library(polyester)
library(Biostrings)

#the command line passed argument is the filepath to the fasta script 
args=commandArgs(trailingOnly=TRUE)

#The polyester library simulates RNAseq reads based on an input fold change matrix and a fasta file containing transcripts to simulate from
#Generate a fold change matrix, with only "1" values since I am just interested in coverage, not differential expression  

primary_fasta_file = args[1]
primary_fasta=readDNAStringSet(primary_fasta_file)

#calculate the number of transcripts to set the number of rows 
sequences=length(primary_fasta)
primary_fold_changes=matrix(c(1), nrow=sequences, ncol=1)

# set the reads_per_transcript: baseline mean number of reads for each transcript
# since I'm just determining mappability, not crucial, but derived this number from the average number of counts per gene and how that related to input reads. 
primary_readspertx=round(20*width(primary_fasta)/150)

#simulate four replicates 
simulate_experiment(primary_fasta_file, reads_per_transcript=primary_readspertx, 
	num_reps=c(2), fold_changes=primary_fold_changes, 
	readlen=75, paired=TRUE, 
	outdir=args[2])
	