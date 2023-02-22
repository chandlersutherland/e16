#download rna reads from ENA 
#using the middle of the 10th leaf 

import os 
import pandas as pd 
import time 

rna_info=pd.read_csv('/global/home/users/chandlersutherland/e16/rna_download_info.txt', sep='\t', header=0)
rna_info['Sample']=rna_info['submitted_ftp'].str.split(pat=';', expand=True)[1].str.split(pat='/', expand=True)[5].str.split(pat='_', expand=True)[0]
rna_info['Tissue']=rna_info['submitted_ftp'].str.split(pat=';', expand=True)[1].str.split(pat='/', expand=True)[5].str.split(pat='_', expand=True)[2]

toi='middle'
subset=rna_info[rna_info['Tissue'] == toi]

for i in range(0, len(subset)):
    start_time = time.time()
    #set sample with for loop
    accession_name=subset.iloc[i,9]
    
    #make dir 
    os.system("mkdir -p /global/scratch/users/chandlersutherland/e16/"+accession_name+"/rna")
    os.chdir("/global/scratch/users/chandlersutherland/e16/"+accession_name+"/rna")
    
    #want the fastq version of the files 
    r1=subset.iloc[i,6].split(';')[0]
    r2=subset.iloc[i,6].split(';')[1]
    
    os.system("wget ftp://"+r1)
    os.system("wget ftp://"+r2)
    print("finished rna download for assembly ", accession_name, ". Total time taken: ", end_time - start_time)
