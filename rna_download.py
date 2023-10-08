#download rna reads from ENA by tissue type 
import os 
import pandas as pd 
import time 
import sys
import glob

#read in metadata file 
rna_info=pd.read_csv('/global/home/users/chandlersutherland/e16/rna_metadata.csv', index_col=0)
sample=sys.argv[1]

#define function that takes the name of a sample, searches to see what has been downloaded, and returns a df of what's left 
def count_left_s(sample):
    #glob paths, convert to df 
    paths=glob.glob(os.path.join('/global/scratch/users/chandlersutherland/e16/'+sample+'/rna_*/*.gz'))
    subset=rna_info[rna_info['Sample']==sample]
    df=pd.DataFrame(paths)
    done=df[0].str.split('/', expand=True)

    #search the subset rna_info file 
    fastq=subset['fastq_ftp'].str.split(';', expand=True)
    subset['fastq_r1']=fastq[0].str.split('/', expand=True)[6]
    subset['fastq_r2']=fastq[1].str.split('/', expand=True)[6]
    
    #melt and search for what's left 
    subset_long=pd.melt(subset, id_vars=['study_accession', 'sample_accession', 'experiment_accession', 
    'run_accession', 'tax_id', 'scientific_name', 'fastq_ftp', 'submitted_ftp', 'sra_ftp', 'Sample', 
    'tissue'],
    value_vars=['fastq_r1', 'fastq_r2'], var_name='fastq_path')
    finished_files=done[8]
    left=subset_long[-subset_long['value'].isin(finished_files)]
    left=left.reset_index()[['study_accession', 'sample_accession', 'experiment_accession', 'run_accession', 'tax_id', 'scientific_name', 'fastq_ftp','submitted_ftp', 'sra_ftp', 'Sample',
    'tissue', 'fastq_path', 'value']]
    to_do=len(left)
    
    #print out 
    print(to_do, 'left in sample ', sample)
    return(left)

#function that takes in that df, then downloads using ftp 
def downloader(left):
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
        if os.path.isfile(file) == False:
            print("download failed for", sample, tissue, file)
        else:
            end_time=time.time()
            print("finished rna download for ", sample, tissue,". Total time taken: ", end_time - start_time)
            

left=count_left_s(sample)
downloader(left)

print('/n/n/n')
print('the following files for sample', sample, 'failed:')
count_left_s(sample) 
