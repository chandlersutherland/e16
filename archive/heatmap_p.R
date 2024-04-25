#goal: create a matrix of p values for hv vs non-hv for different RNA tissue types and methylation states 
library(pheatmap)

p_values <- read_csv('//wsl.localhost//Ubuntu//home//chandlersutherland//e16_scratch//p_value_matrix.csv', )
p_values$accession <- subpopulations$accession %>% toupper()
before_correction <- column_to_rownames(p_values, 'accession') %>% subset(select=c('CpG', 'CHH', 'CHG'))
after_correction <- column_to_rownames(p_values, 'accession') %>% subset(select=c('CpG_BH', 'CHH_BH', 'CHG_BH'))
before_correction
myBreaks <- c(0, 0.001, 0.01, 0.05, 1)
myCols <- c('red', 'orange', 'yellow', 'grey')
pheatmap(before_correction, breaks=myBreaks, color=myCols)
pheatmap(after_correction, breaks=myBreaks, color=myCols, cluster_cols=FALSE)

#methylation only, correct BH across rows not columns (independent measurements)
meth <- before_correction %>% subset(select=c('CpG', 'CHH', 'CHG'))
p.adjust(meth[1,], method='BH')
new <- apply(meth, 1, p.adjust) %>% t()
pheatmap(new, breaks=myBreaks, color=myCols, cluster_cols=FALSE)

#subpopulation info as an annotation column 
subpopulations <- read_table('//wsl.localhost//Ubuntu//home//chandlersutherland//e16//nam_genome_info.txt', col_names=c('Assembly', 'Grin', 'accession_id', 'source', 'cross_reference', 'subpopulation', 'stock'), skip=1) %>% separate(Assembly, sep='-', into=c(NA, 'accession', NA, NA, NA))
subpopulations$accession <- subpopulations$accession %>% toupper()
subpopulations <- subpopulations %>% mutate(subpopulation =recode(subpopulation, 'Temporate/tropical'='Temporate/Tropical'))%>% 
  mutate(subpopulation =recode(subpopulation, 'Temporate/Tropical'='Mixed', 'Non-stiff-stalk'='Non-stiff stalk', 'Sweet'='Sweetcorn'))
df <- subpopulations %>% subset(select=c('accession', 'subpopulation')) %>% column_to_rownames('accession')
df$subpopulation <- factor(df$subpopulation, levels=c('Stiff stalk', 'Non-stiff stalk', 'Mixed', 'Popcorn', 'Sweetcorn', 'Tropical'))
my_col = list(
  subpopulation = c('Non-stiff stalk'='#80d6f6', 
                      'Mixed'='#d7d7d7', 
                      'Popcorn'='#e6b1ff', 
                      'Sweetcorn'='#ffbf70', 
                      'Tropical'='#78d5a0'))

met1 <- pheatmap(new, breaks=myBreaks, color=myCols, cluster_cols=FALSE, annotation_row=df, annotation_colors=my_col)
ggsave('C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//met_p.png', plot=met1, dpi='retina', width=5, height=5)

tip <- p_values %>% column_to_rownames('accession') %>% subset(select=c('tip'))
tip_p <- pheatmap(tip, breaks=myBreaks, color=myCols, cluster_cols=FALSE, annotation_row=df, annotation_colors=my_col)
ggsave('C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//tip_p.png', plot=tip_p, dpi='retina', width=5, height=5)


#repeat with all tissue info as of 10.17
morep_values <- read_csv('//wsl.localhost//Ubuntu//home//chandlersutherland//e16_scratch//tissue_p_matrix.csv', )
morep_values$accession <- morep_values$accession %>% toupper()
before_correction <- column_to_rownames(morep_values, 'accession')
morep_values
pheatmap(before_correction, breaks=myBreaks, color=myCols)
#l<-sapply(before_correction, function(x) p.adjust(x, method="BH"))
#after_correction <- column_to_rownames(p_values, 'accession') %>% lapply(p.adjust(method='BH'))

morepbh_values <- read_csv('//wsl.localhost//Ubuntu//home//chandlersutherland//e16_scratch//tissue_pbh_matrix.csv', )
morepbh_values$accession <- morepbh_values$accession %>% toupper()
after_correction <- column_to_rownames(morepbh_values, 'accession')

bh <- pheatmap(after_correction, breaks=myBreaks, color=myCols, annotation_row=df, annotation_colors=my_col)
ggsave('C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//tissue_p.png', plot=bh, dpi='retina', width=7, height=5)

