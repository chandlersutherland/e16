##---
##title: "all_hv_enrich"
##author: "Chandler Sutherland"
##date: "2023-11-15"
##---
setwd('/global/scratch/users/chandlersutherland/e16/cs_reports')

library(tidyverse)
library(ggplot2)
library(data.table)
library(singscore)
library(GSEABase)
library(patchwork)
library(pheatmap)
#library(stringi)

#step 1: read in all hv file, clean up 
all_hv <- read_lines('Zm_hv_OG_geneIDs.txt')
#test <- all_hv[1] %>% str_split(' ')
#clade <- test[[1]][1] %>% str_replace(':', '')
#length(test[[1]])
#genes <- test[[1]][2:8]
#genes
#df <- tibble(genes) %>% mutate(clade=clade)

#length(all_hv)

converter <- function(line){
  test <- line %>% str_split(' ')
  all_hv_clade <- test[[1]][1] %>% str_replace(':', '')
  l <- length(test[[1]]) -1 
  Gene <- test[[1]][2:l]
  df <- tibble(Gene) %>% mutate(Clade=all_hv_clade)
  return(df) 
}

res <- lapply(all_hv, converter)
all_hv_df <- as.data.frame(do.call(rbind, res)) %>% separate(Gene, '_', into=c("Gene", NA))
all_hv_genes <- all_hv_df$Gene 

#step 2: read in all tpm files 
tpm_files <- Sys.glob('/global/scratch/users/chandlersutherland/e16/cs_reports/*_all_tissue.tsv')

read_in <- function(file_path){
  f <- read.csv(file_path, sep='\t') %>% subset(select=c('accession', 'tissue',  'rep', 'name', 'TPM', 'stranded_1'))
  return(f)
}

#step 2.1: filter lowcount genes 
filter_low <- function(f, min_count, percent_samples){
  #min_count is min TPM across percent samles
  expr <- f %>% group_by(name) %>% summarize(top=quantile(TPM, percent_samples))%>%
    filter(top > min_count) %>% pull(name)
  filtered <- f %>% filter(name %in% expr)  
  
  return(filtered)
}



b73 <- read_in(tpm_files[1])
b73_filter <- filter_low(b73, 10, .9)
split_b73 <- split(b73_filter, b73_filter$tissue)
test <- split_b73[[1]]

ggplot(b73, aes(x=log2(TPM+1)))+
  geom_histogram(aes(color=tissue))

expr <- b73 %>% group_by(name) %>% summarize(top=quantile(TPM, .9))%>%
  filter(top > 5) %>% pull(name)
b73_filter <- b73 %>% filter(name %in% expr) 

ggplot(b73_filter, aes(x=log2(TPM+1)))+
  geom_histogram(aes(color=tissue))
#step 3: create hv gene set 
hv_geneset_id <- function(split_df){
  hv_names <- split_df$name[split_df$name %in% all_hv_genes] %>% unique()
  hv_geneset <- GSEABase::GeneSet(hv_names, geneIdType=SymbolIdentifier())
  return(hv_geneset)
}

#step 4: build tpm matrix 
ranked <- function(split_df){
  matrix <- split_df %>% 
    subset(select=c('rep', 'name', 'TPM')) %>% 
    pivot_wider(names_from=rep, values_from=TPM)
  matrix <- column_to_rownames(matrix, 'name')
  ranked <- rankGenes(matrix)
  return(ranked)
}

#step 5: calculate p value of hv set 
singscore_p <- function(ranked, hv_geneset){
  hv_sing <- simpleScore(ranked, upSet=hv_geneset)
  hv_null <- generateNull(rankData=ranked, upSet=hv_geneset)
  hv_p <- getPvals(hv_null, hv_sing) #%>% as.data.frame()
  #res <- as.data.frame(hv_p) %>% 
   # rownames_to_column('rep') %>% 
    #mutate(tissue=tissue) %>% 
    #mutate(accession=accession)
  return(hv_p)
}

l <- hv_geneset_id(test)
p <- ranked(test)
res <- singscore_p(p, l)

#step 6: plot 
singscore_plot <- function(ranked, hv_geneset, tissue){
  n_plot <- ncol(ranked)[[1]]
  reps <- colnames(ranked)
  plots <- plot_spacer()
  for (x in 1:n_plot) {
    rep <- colnames(ranked)[x]
    p <- plotRankDensity(ranked[,x,drop=FALSE], upSet=hv_geneset, isInteractive=FALSE)+
      xlab('')+
      ylab('')+
      ggtitle(paste(tissue,reps[x]))+
      theme(text=element_text(size=5))
    plots <- plots | p + theme(text=element_text(size=5))
  }
  #plots <- plots + plot_annotation(tissue)
  return(plots)
}
#name(ranked(test))
t <- singscore_plot(ranked(test), hv_geneset_id(test), 'test')

split_b73
all_sets <- lapply(split_b73, hv_geneset_id)
all_ranked <- lapply(split_b73, ranked)
r <- mapply(singscore_p, ranked=all_ranked, hv_geneset=all_sets)
r_df <- r %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rownames_to_column('tissue') %>% 
  separate(tissue, into=c('tissue', 'rep')) 

colnames(r_df) <- c('tissue', 'rep', 'singscore_p')
p <- mapply(singscore_plot, ranked=all_ranked, hv_geneset=all_sets, tissue=names(all_ranked))
t <- wrap_plots(p, ncol=2) + theme(text=element_text(size=5)) + plot_annotation('B73')
filepath=paste('/global/scratch/users/chandlersutherland/e16/cs_reports/hv_enrichment_plots/', 'B73', '.png', sep='')
ggsave(filepath, t, width=11, height=8.5)

#step 7: wrap 
wrapper <- function(in_file_path){
  #get accession name from file path  
  b <- in_file_path %>% str_split('/')
  f_name <- b[[1]][8] %>% str_split('_')
  accession <- f_name[[1]][1]
  print(paste('starting accession:', accession))
  
  #input and filter data
  f <- read_in(in_file_path)
  filtered <- filter_low(f, 5, .9)
  split_df <- split(filtered, filtered$tissue)
  
  #run singscore 
  all_sets <- lapply(split_df, hv_geneset_id)
  all_ranked <- lapply(split_df, ranked)
  r <- mapply(singscore_p, ranked=all_ranked, hv_geneset=all_sets)
  
  #add processing of r 
  r_df <- r %>% 
    unlist() %>% 
    as.data.frame() %>% 
    rownames_to_column('tissue') %>% 
    separate(tissue, into=c('tissue', 'rep')) 
  colnames(r_df) <- c('tissue', 'rep', 'singscore_p')
  r_df$accession <- accession 
  
  #plot and save to file 
  p <- mapply(singscore_plot, ranked=all_ranked, hv_geneset=all_sets, tissue=names(all_ranked))
  t <- wrap_plots(p, ncol=2) + 
    theme(text=element_text(size=5)) + 
    plot_annotation(accession)
  filepath=paste('/global/scratch/users/chandlersutherland/e16/cs_reports/hv_enrichment_plots/', accession, '.png', sep='')
  ggsave(filepath, t, width=11, height=8.5)
  
  #return p value df 
  return(r_df)
}

wrapper(tpm_files[1])
all_p <- lapply(tpm_files, wrapper) #this takes ~4 minutes per accession, 2hr total

#create p value df, and write to csv 

#bh correct p values 

#create a matrix with accessions as rows and tissue type as column, largest p value of each sample as value 



#plot heatmap after bh correction 
myBreaks <- c(0, 0.001, 0.01, 0.05, 1)
myCols <- c('red', 'orange', 'yellow', 'grey')
pheatmap(all_p, breaks=myBreaks, color=myCols)

p39 <- read_in(tpm_files[24])
filtered <- filter_low(p39, 5, .9)
split_df <- split(filtered, filtered$tissue)
all_sets <- lapply(split_df, hv_geneset_id)
all_ranked <- lapply(split_df, ranked)
r <- mapply(singscore_p, ranked=all_ranked, hv_geneset=all_sets)
#if you have a perfect two samples per replicate, the current treatment of r df doesn't work 

