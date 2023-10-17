#download rna reads from ENA by tissue type 
import os 
import pandas as pd 
import time 
import sys
import glob

#read in metadata file 
rna_info=pd.read_csv('/global/home/users/chandlersutherland/e16/rna_metadata.csv', index_col=0)

#supply accession as sysargv[1]
#sample=sys.argv[1]
#subset=rna_info[rna_info['Sample']==sample]

#for i in range(0, len(subset)):
 #   start_time = time.time()
    
    #set tissue 
#    tissue=subset.iloc[i,10]
    
    #make dir 
#    os.system("mkdir -p /global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
#    os.chdir("/global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
    
    #want the fastq version of the files 
#    r1=subset.iloc[i,6].split(';')[0]
#    r2=subset.iloc[i,6].split(';')[1]
    
#    os.system("wget ftp://"+r1)
#    os.system("wget ftp://"+r2)
#    end_time=time.time()
#    print("finished rna download for ", sample, tissue,". Total time taken: ", end_time - start_time)
    
#restarting, check what is there and ask to find new stuff again... 
#system argument should be tissue type 
tissue=sys.argv[1]
#print(tissue)
#collect paths of existing downloads
paths=glob.glob(os.path.join('/global/scratch/users/chandlersutherland/e16/*/rna_'+tissue+'/*.gz'))
#print(paths)
subset=rna_info[rna_info['tissue']==tissue]

#make into a dataframe 
df=pd.DataFrame(paths)
done=df[0].str.split('/', expand=True)

#make rna_info subset longer to account for both rna read files 
fastq=subset['fastq_ftp'].str.split(';', expand=True)
subset['fastq_r1']=fastq[0].str.split('/', expand=True)[6]
subset['fastq_r2']=fastq[1].str.split('/', expand=True)[6]

#melt df 
subset_long=pd.melt(subset, id_vars=['study_accession', 'sample_accession', 'experiment_accession', 'run_accession', 'tax_id', 'scientific_name', 'fastq_ftp', 'submitted_ftp', 'sra_ftp', 'Sample', 'tissue'],
 value_vars=['fastq_r1', 'fastq_r2'], var_name='fastq_path')

#get just remaining files to download     
finished_files=done[8]
left=subset_long[-subset_long['value'].isin(finished_files)]
left=left.reset_index()[['study_accession', 'sample_accession', 'experiment_accession', 'run_accession', 'tax_id', 'scientific_name', 'fastq_ftp','submitted_ftp', 'sra_ftp', 'Sample', 'tissue', 'fastq_path', 'value']]

len_left=len(left)
print(str(len_left)+'files left to download for'+tissue)

#finally restart for loop 
for i in range(0, len(left)):
    start_time = time.time()
    sample=left.iloc[i, 9]
    tissue=left.iloc[i, 10]
    os.system("mkdir -p /global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
    os.chdir("/global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
    
    r1=left.iloc[i,6].split(';')[0]
    r2=left.iloc[i,6].split(';')[1]
    
    os.system("wget ftp://"+r1)
    os.system("wget ftp://"+r2)
    
    file=left.iloc[i,12]
    if os.path.isfile(file) == True:
        print("download failed for", sample, tissue, file)
    else: 
        end_time=time.time()
        print("finished rna download for ", sample, tissue,". Total time taken: ", end_time - start_time)
    
