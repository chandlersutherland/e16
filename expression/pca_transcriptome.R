library(data.table)
library(ggplot2)
library(tidyverse)
library(purrr)
library(viridis)

#title: PCA of variance transcriptome wide 
#author: chandler sutherland

#step 1: read in pangene matrix that was released with hufford et al 
pan_gene <- read.csv('/global/scratch/users/chandlersutherland/e16/cs_reports/pan_gene_matrix_v3_cyverse.csv')
core <- pan_gene %>% filter(class=='Core Gene') #27k core genes 
uppercase <- colnames(core)[4:29] %>% toupper()
new_names <- append(uppercase, c('PAN_GENE_ID', 'UNIVERSAL_GENE_ID', 'Rep_transcript'), after=0) %>% append(c('class', 'Subgenome', 'number_genome_presence'))
colnames(core) <- new_names 

#step 2: read in tpm files for all tissues 
tpm_files <- Sys.glob('/global/scratch/users/chandlersutherland/e16/cs_reports/*_all_tissue.tsv')

#define a function that reads in, calculates mean log2TPM, and returns a tidy df for each accession 
read_in <- function(file_path){
  #import tsv file 
  f <- read.csv(file_path, sep='\t')
  
  #calculate mean TPM, not interested in between rep variation we know they cluster together
  f_reformat <- f %>% group_by(accession, tissue, name) %>% 
    summarize(log2_TPM=mean(log2.TPM., na.rm=TRUE)) %>% 
                subset(select=c('accession','tissue',  'name', 'log2_TPM'))
  return(f_reformat)
}

#define a function that converts the df from read_in into a matrix with columns with pan gene id, and rows with accession+tissue
matrix_maker <- function(f_reformat){
  acc <- f_reformat[[1,1]]
  print(acc)
  c <- core %>% subset(select=c('PAN_GENE_ID', acc)) %>% separate(acc, sep='_', into=c(acc, NA)) #subset and clean up core 
  named <- merge(c, f_reformat, by.x=acc, by.y='name') #add pan gene IDs
  m <- named %>% mutate(col_id=paste(accession, tissue, sep='_')) %>% 
    subset(select=c('PAN_GENE_ID', 'col_id', 'log2_TPM')) %>% 
    pivot_wider(names_from=PAN_GENE_ID, values_from=log2_TPM)
  return(m)
}

#perform read_in 
raw <- lapply(tpm_files, read_in)

#perform matrix_maker
matrices <- lapply(raw, matrix_maker)

#create the master matrix 
master_matrix <- rbindlist(matrices, use.names=T, fill=T)
master_matrix %>% write.csv('/global/scratch/users/chandlersutherland/e16/cs_reports/master_tpm_matrix.csv')

#read in saved copy to save compute time 
master_matrix <- read_csv('/global/scratch/users/chandlersutherland/e16/cs_reports/master_tpm_matrix.csv') %>% subset(select=-c(`...1`))

#QC matrix 
ncol(master_matrix) #recover 27,910 pan genes 
na_count <- master_matrix %>% summarise(across(everything(), ~ sum(is.na(.))))
tidy <- na_count %>% t() %>% as.data.frame()
colnames(tidy) <- 'count'
ggplot(tidy)+geom_histogram(aes(x=count))
tidy %>% filter(count > 0) %>% nrow() #6,472 genes have someone missing.. 
missing <- tidy %>% filter(count > 0) %>% rownames_to_column('PAN_GENE_ID') %>% pull('PAN_GENE_ID')
found <- tidy %>% filter(count == 0) %>% rownames_to_column('PAN_GENE_ID') %>% pull('PAN_GENE_ID')
core %>% filter(PAN_GENE_ID %in% missing)
#missing genes are coming from these funky gmap coded genes. 
#find TZI8 gmap_ID=chr1:1117701-1120773
TZI8 <- read.csv('/global/scratch/users/chandlersutherland/e16/cs_reports/TZI8_all_tissue.tsv', sep='\t') 
TZI8 %>% filter(chrom=='chr1') %>% filter(chromStart == 1117701)

#calculate mean TPM, not interested in between rep variation we know they cluster together
f_reformat <- f %>% group_by(accession, tissue, name) %>% 
  summarize(log2_TPM=mean(log2.TPM., na.rm=TRUE))

#create PCA plot! 
pca_matrix <- master_matrix %>%
  subset(select=found)%>% #filter any pangene with na values 
  column_to_rownames('col_id') %>% 
  as.matrix() #%>% 
  #na.omit()

sample_pca <- prcomp(pca_matrix)
head(sample_pca)
#examine variance explained by PCAs 
pc_eigenvalues <- sample_pca$sdev^2 #vector with variance explained by each PC and how much variance explained 

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#get PC scores 
pc_scores <- sample_pca$x
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")  %>% 
  separate(sample, into=c('accession', 'tissue'))

# print the result
pc_scores %>% write.csv('/global/scratch/users/chandlersutherland/e16/cs_reports/PC_tpm_matrix.csv')

#check what the result would be without endosperm/embryo, which have missing representative samples
expr <- master_matrix
subsample <- master_matrix %>% filter(!str_detect(col_id, "_endosperm")) %>% filter(!str_detect(col_id, '_embryo'))
colSums(subsample, na.rm=TRUE)
pca_matrix2 <- subsample %>%
  subset(select=found)%>% #filter any pangene with na values 
  column_to_rownames('col_id') #%>% 
  #as.matrix()
sum(colSums(pca_matrix2, na.rm=TRUE) == 0)
sample_pca2 <- prcomp(pca_matrix2)
head(sample_pca2)
#examine variance explained by PCAs 
pc_eigenvalues2 <- sample_pca2$sdev^2 #vector with variance explained by each PC and how much variance explained 

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues2 <- tibble(PC = factor(1:length(pc_eigenvalues2)), 
                         variance = pc_eigenvalues2) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues2

pc_eigenvalues2 %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#get PC scores 
pc_scores2 <- sample_pca2$x
pc_scores2 <- pc_scores2 %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")  %>% 
  separate(sample, into=c('accession', 'tissue'))

pc_scores2 %>% write.csv('/global/scratch/users/chandlersutherland/e16/cs_reports/PC2_tpm_matrix.csv')

pc_scores2 %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, color=tissue)) +
  geom_point(size=0.5)+
  xlab('PC1: 20.2% of variance')+
  ylab('PC2: 10.9% of variance')+
  theme_classic()+
  scale_color_viridis(discrete=TRUE, option="turbo")+
  labs(color='Tissue')+
  theme(text=element_text(size=10))+
  # theme(legend.spacing.y = unit(.05, 'cm')) +
  theme(legend.key.height=unit(.8,"line"))+
  guides(fill = guide_legend(byrow = TRUE))

