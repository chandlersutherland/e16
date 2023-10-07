#download rna reads from ENA by tissue type 
import os 
import pandas as pd 
import time 
import sys

#read in metadata file 
rna_info=pd.read_csv('/global/home/users/chandlersutherland/e16/rna_metadata.csv', index_col=0)

#supply accession as sysargv[1]
sample=sys.argv[1]
subset=rna_info[rna_info['Sample']==sample]

for i in range(0, len(subset)):
    start_time = time.time()
    
    #set tissue 
    tissue=subset.iloc[i,10]
    
    #make dir 
    os.system("mkdir -p /global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
    os.chdir("/global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
    
    #want the fastq version of the files 
    r1=subset.iloc[i,6].split(';')[0]
    r2=subset.iloc[i,6].split(';')[1]
    
    os.system("wget ftp://"+r1)
    os.system("wget ftp://"+r2)
    end_time=time.time()
    print("finished rna download for ", sample, tissue,". Total time taken: ", end_time - start_time)