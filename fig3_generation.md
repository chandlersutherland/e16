fig_generation
================
Chandler Sutherland
2023-11-17

Purpose: Generating the plots shown in Fig3 and Supplemental Figure 3

``` r
library(data.table)
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.3.1

``` r
library(tidyverse)
```

    ## Warning: package 'purrr' was built under R version 4.3.1

    ## Warning: package 'dplyr' was built under R version 4.3.1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.0
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::between()     masks data.table::between()
    ## ✖ dplyr::filter()      masks stats::filter()
    ## ✖ dplyr::first()       masks data.table::first()
    ## ✖ lubridate::hour()    masks data.table::hour()
    ## ✖ lubridate::isoweek() masks data.table::isoweek()
    ## ✖ dplyr::lag()         masks stats::lag()
    ## ✖ dplyr::last()        masks data.table::last()
    ## ✖ lubridate::mday()    masks data.table::mday()
    ## ✖ lubridate::minute()  masks data.table::minute()
    ## ✖ lubridate::month()   masks data.table::month()
    ## ✖ lubridate::quarter() masks data.table::quarter()
    ## ✖ lubridate::second()  masks data.table::second()
    ## ✖ purrr::transpose()   masks data.table::transpose()
    ## ✖ lubridate::wday()    masks data.table::wday()
    ## ✖ lubridate::week()    masks data.table::week()
    ## ✖ lubridate::yday()    masks data.table::yday()
    ## ✖ lubridate::year()    masks data.table::year()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggsignif)
library(ggpubr)
library(ggbeeswarm)
library(pheatmap)
```

    ## Warning: package 'pheatmap' was built under R version 4.3.1

``` r
library(patchwork)
```

    ## Warning: package 'patchwork' was built under R version 4.3.1

``` r
library(introdataviz)
library(viridis)
```

    ## Warning: package 'viridis' was built under R version 4.3.1

    ## Loading required package: viridisLite

Import and clean data

``` r
#remove unmappable NLRs 
unmappable <- c('Zm00039ab351270', 'Zm00026ab135540', 'Zm00036ab418650',
                'Zm00001eb091500', 'Zm00001eb164880', 'Zm00001eb343890',
                'Zm00001eb391100', 'Zm00001eb405870', 'Zm00001eb405900',
                'Zm00001eb405930', 'Zm00033ab429000', 'Zm00029ab367660')

#import gene information 
gene_table <- read.csv('Maize_NLRome_GeneTable.txt', sep='\t') %>% subset(select=c('Gene', 'Ecotype', 'HV', 'Clade'))
gene_table$Ecotype <- gene_table$Ecotype %>% toupper()
gene_table$Gene <- gene_table$Gene %>% 
  str_replace('ZM', 'Zm') %>% 
  str_replace('AB', 'ab') %>% 
  str_replace("EB", 'eb') %>% 
  str_replace('_P001', '')

#import nlr TPM table, average biological replicates 
all_tpm <- read.csv('//wsl.localhost//Ubuntu//home//chandlersutherland//e16_scratch//all_nlr_tissue.csv') %>% 
  filter(!tissue %in% c('endosperm', 'embryo')) %>% 
  filter(!name %in% unmappable) %>% 
  merge(gene_table, by.x='name', by.y='Gene')

all_tpm_avg <- all_tpm %>%
  group_by(Clade, HV, accession, tissue, name, chrom, chromStart, chromEnd, strand) %>%
  summarize(log2_TPM=mean(log2.TPM., na.rm = T)) %>%
  ungroup() 
```

    ## `summarise()` has grouped output by 'Clade', 'HV', 'accession', 'tissue',
    ## 'name', 'chrom', 'chromStart', 'chromEnd'. You can override using the `.groups`
    ## argument.

``` r
#calculate significance of difference, correct by Benjamani-Hochberg 
p_tbl <- compare_means(log2_TPM~HV, data=all_tpm_avg, group.by='tissue', p.adjust.method = 'BH') %>% arrange(p.adj)
order_p <- p_tbl %>% pull(tissue)
p_annot <- p_tbl %>% pull(p.signif)

#clean up data for pretty plotting 
all_tpm_avg$tissue <- factor(all_tpm_avg$tissue, levels=c(rev(order_p)))
all_tpm_avg <- all_tpm_avg %>% mutate(HV=case_match(HV, 0 ~ "non-hv", 1~"hv")) #I don't know why this takes so long.. 
all_tpm_avg$HV <- factor(all_tpm_avg$HV, levels=c('non-hv', 'hv'))
```

Create an expression matrix

``` r
#reformat tpm table 
#define a function that applies a sample ID based on the replicates
reformat <- function(split_df){
  #by_tissue <- split(split_df, split_df$tissue)
  split_df <- split_df  %>% group_by(rep) %>% mutate(sample_id=cur_group_id()) %>% mutate(tissue=paste(tissue, as.character(sample_id), sep='_'))
  split_df 
}

#define a function that applies reformat to each accession data frame 
by_tissue <- function(by_accession){
  one_acc <- lapply(split(by_accession, by_accession$tissue), reformat) %>% bind_rows()
  one_acc
}

#bring it all together 
all <- all_tpm  %>%  split(all_tpm$accession)
reformed <- lapply(all, by_tissue) %>% bind_rows() %>% subset(select=c(HV, accession, name, tissue, log2.TPM.)) %>%
  unique() %>%
  pivot_wider(names_from=tissue, values_from=log2.TPM.) %>% 
  remove_rownames() %>% 
  column_to_rownames('name')

#create an hv expression matrix 
p <- reformed %>% 
  filter(HV==1) %>% 
  subset(select=-c(HV, accession, tip_3, middle_3, root_3, base_3)) %>% #remove extra bioreps 
  filter_all(any_vars(. != 0))  %>% #remove genes with no variation across tissue types (would also remove all 0 genes)
  subset(select=c('root_1', 'root_2', 'ear_1', 'ear_2', 'shoot_1', 'shoot_2', 'base_1', 'base_2', 'tassel_1', 'tassel_2',
                  'anther_1', 'anther_2', 
                   'middle_1', 'middle_2', 'tip_1', 'tip_2'))
```

To define the gene expression categories, we compared TPM values across
the 10 tissues (leaf tip, leaf middle, leaf base, root, shoot, ear,
anther, tassel, endosperm, and embryo). We defined tissue-specific
expression as TPM ≥ 1 in at least 1 tissue and TPM \< 1 in at least 1
tissue, constitutive expression as TPM ≥ 1 in all 10 tissues, and silent
as TPM \< 1 in all 10 tissues.

``` r
by_category <- all_tpm %>% 
  group_by(accession, name, tissue, HV, chrom, chromStart,chromEnd, strand) %>%
  summarize(TPM=mean(TPM, na.rm=TRUE)) %>%
  mutate(HV=case_match(HV, 0 ~ "non-hv", 1~"hv")) %>%
  pivot_wider(names_from=tissue, values_from=TPM) %>%
  mutate(expr_category=case_when((
    (
      (anther >= 1) | (base >= 1) | (ear >= 1) | (middle >= 1) | (root >= 1) | (shoot >= 1) | (tip >= 1) | (tassel >=1)
      ) & 
      (
        (anther < 1) | (base < 1) | (ear < 1) | (middle < 1) | (root < 1) | (shoot < 1) | (tip < 1) | (tassel < 1)
        )
    )  ~'tissue_specific',
    (
      (anther >= 1) & (base >= 1) & (ear >= 1) & (middle >= 1) & (root >= 1) & (shoot >= 1) & (tip >= 1) & (tassel >=1)
      ) ~ 'constitutive',
    (
      (anther < 1) & (base < 1) & (ear < 1) & (middle < 1) & (root < 1) & (shoot < 1) & (tip < 1) & (tassel < 1)
      ) ~ 'silent'))
```

    ## `summarise()` has grouped output by 'accession', 'name', 'tissue', 'HV',
    ## 'chrom', 'chromStart', 'chromEnd'. You can override using the `.groups`
    ## argument.

``` r
by_category$HV <- factor(by_category$HV, levels=c('non-hv','hv'))


all_tpm_avg <- by_category %>% subset(select=c(name, expr_category)) %>% unique() %>% merge(all_tpm_avg) %>% 
  mutate(Clade_adj=recode(Clade, 'Int3480_75_130_L_68'='Int3480', 
                          'Int4787_129_172_L_83'='RppM-like', 
                          'Int6329_131_253_L_122_209_L_80_142_L_61'='RppC-like', 'Int6329_131_253_L_122_209_R_35_43_R_23'='RppC-like',
                          'Int6329_131_253_L_122_209_L_80_142_R_18'='RppC-like',
                          'Int6648_150_178_L_128_144_R_93'='Rp1-like'))

write_csv(all_tpm_avg, file="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//intermediate_data//nlr_tpm.csv")

coordinates <- all_tpm_avg %>% 
  subset(select=c('accession', 'name', 'Clade_adj', 'HV', 'expr_category', 'chrom', 'chromStart', 'chromEnd', 'strand')) %>% 
  unique()

write_csv(coordinates, file="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//intermediate_data//nlr_coords.csv")
```

## Fig 3A: per clade percentage category

``` r
all_tpm_avg$expr_category <- factor(all_tpm_avg$expr_category, levels=c('silent', 'tissue_specific', 'constitutive'))  

cat_plot <- all_tpm_avg %>% 
  group_by(HV, Clade_adj, expr_category) %>% 
  summarize(n=n()) %>% 
  mutate(prop=n/sum(n)*100) %>% 
  ungroup()%>%
  group_by(HV, expr_category) %>% 
  mutate(rank=rank(prop)) %>%
  mutate(Clade_adj=reorder(Clade_adj, -rank)) %>%
  ungroup()
```

    ## `summarise()` has grouped output by 'HV', 'Clade_adj'. You can override using
    ## the `.groups` argument.

``` r
#sample size of clades 
sample_size <- cat_plot %>% subset(select=c(HV, Clade_adj)) %>% unique() %>% group_by(HV) %>% summarize(n=n()) %>% pull(n)
xlabel <- paste(c('Clade\nn=',
                  as.character(sample_size[[1]]), 
                  ' non-hvNLRs, ', 
                  as.character(sample_size[[2]]), 
                  ' hvNLRs'), 
                sep='', collapse='')

#print order of hv clades 
cat_plot %>% filter(HV=='hv') %>% pull(Clade_adj) %>% unique()
```

    ## [1] Int3480   Rp1-like  RppC-like RppM-like
    ## 188 Levels: Int3452_26 Int3480_75_130_R_1 ... RppC-like

``` r
by_category_plot <- cat_plot  %>%
  ggplot(aes(x=Clade_adj, y=prop, fill=expr_category))+
  geom_col()+
  facet_grid(~HV, scales='free', space='free')+
  ylab('% of Clade')+
  xlab(xlabel)+
  #labs(fill='Gene Expression Category')+
  scale_fill_discrete(name = "Gene Expression Category", labels = c("Silent", "Tissue Specific", "Constitutive"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=10)
        #axis.ticks.x = element_blank()
        )


by_category_plot 
```

![](fig3_generation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave(by_category_plot, filename='C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//expression_by_category.svg', width=180, height=70, units='mm')
```

## Fig 3B: Heatmap

``` r
#create an annotation matrix 

#import subpopulation membership 
subpopulations <- read_table('//wsl.localhost//Ubuntu//home//chandlersutherland//e16//nam_genome_info.txt', col_names=c('Assembly', 'Grin', 'accession_id', 'source', 'cross_reference', 'subpopulation', 'stock'), skip=1) %>% separate(Assembly, sep='-', into=c(NA, 'accession', NA, NA, NA))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Assembly = col_character(),
    ##   Grin = col_character(),
    ##   accession_id = col_double(),
    ##   source = col_character(),
    ##   cross_reference = col_character(),
    ##   subpopulation = col_character(),
    ##   stock = col_character()
    ## )

    ## Warning: 25 parsing failures.
    ## row col  expected    actual                                                                          file
    ##   1  -- 7 columns 8 columns '//wsl.localhost//Ubuntu//home//chandlersutherland//e16//nam_genome_info.txt'
    ##   2  -- 7 columns 8 columns '//wsl.localhost//Ubuntu//home//chandlersutherland//e16//nam_genome_info.txt'
    ##   3  -- 7 columns 8 columns '//wsl.localhost//Ubuntu//home//chandlersutherland//e16//nam_genome_info.txt'
    ##   4  -- 7 columns 8 columns '//wsl.localhost//Ubuntu//home//chandlersutherland//e16//nam_genome_info.txt'
    ##   5  -- 7 columns 8 columns '//wsl.localhost//Ubuntu//home//chandlersutherland//e16//nam_genome_info.txt'
    ## ... ... ......... ......... .............................................................................
    ## See problems(...) for more details.

``` r
subpopulations$accession <- subpopulations$accession %>% toupper()
subpopulations <- subpopulations %>% subset(select=c('accession', 'subpopulation'))
subpopulations[nrow(subpopulations) + 1,] = list('B73', 'Stiff stalk')
subpopulations <- subpopulations %>% mutate(subpopulation =recode(subpopulation, 'Temporate/tropical'='Temporate/Tropical'))%>% 
  mutate(subpopulation =recode(subpopulation, 'Temporate/Tropical'='Mixed', 'Non-stiff-stalk'='Non-stiff stalk', 'Sweet'='Sweetcorn'))
subpopulations$subpopulation <- factor(subpopulations$subpopulation, levels=c('Stiff stalk', 'Non-stiff stalk', 'Mixed', 'Popcorn', 'Sweetcorn', 'Tropical'))


#apply clade membership from geneTable
row_annotation <- all_tpm_avg %>%
  filter(HV=='hv') %>% 
  subset(select=c(name, accession, Clade_adj, expr_category)) %>% 
  unique() %>% 
  merge(subpopulations) %>% 
  subset(select=c(name,  
                  expr_category, 
                  #subpopulation, 
                  Clade_adj)) %>% 
  remove_rownames() %>% 
  column_to_rownames('name')

#create column annotations
samples <- colnames(p)
col_annotation <- as.data.frame(samples) %>% separate(samples, '_', into=c('tissue', NA), remove=F) %>% column_to_rownames('samples')

#create color annotation list 
my_col = list(
  Clade_adj=c('Int3480'='#fee090', 
                              'RppM-like'='#d73027', 
                              'RppC-like'='#fc8d59', 
                              'Rp1-like'='#4575b4'), 
  #subpopulation=c('Stiff stalk'='#ffe3a9', 
  #                            'Non-stiff stalk'='#80d6f6', 
  #                            'Mixed'='#d7d7d7', 
  #                            'Popcorn'='#e6b1ff', 
  #                            'Sweetcorn'='#ffbf70', 
  #                            'Tropical'='#78d5a0'),
  expr_category=c('tissue_specific'='#00BA38', 
                  'silent'= '#F8766D',
                  'constitutive'='#619CFF'
                  ),
  tissue=c('anther'='#30123BFF', 
           'base'='#4777EFFF', 
           'ear'='#1BD0D5FF', 
           'middle'='#62FC6BFF', 
           'root'='#D2E935FF', 
           'shoot'='#FE9B2DFF', 
           'tassel'='#DB3A07FF', 
           'tip'='#7A0403FF'))

#get gene lists for each hv clade
Int3480 <- row_annotation %>% filter(Clade_adj=='Int3480') %>% rownames()
Rp1 <- row_annotation %>% filter(Clade_adj=='Rp1-like') %>% rownames()
RppM <- row_annotation %>% filter(Clade_adj=='RppM-like') %>% rownames()
RppC <- row_annotation %>% filter(Clade_adj=='RppC-like') %>% rownames()

#standardize mat breaks to I can create each heatmap individually 
mat_breaks <- seq(min(p, na.rm=TRUE), max(p, na.rm=TRUE), length.out = 100)
```

``` r
#plot! 
Int3480_m <- p[rownames(p) %in% Int3480,]
pheatmap(Int3480_m, 
         cluster_cols=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE,
         annotation_row=row_annotation, 
         annotation_colors=my_col, 
         cellheight=1.25, cellwidth=5, 
         border_color=NA,
         legend=FALSE,
         annotation_legend=FALSE, 
         annotation_names_row=FALSE,
         treeheight_row=30,
         breaks=mat_breaks, 
        filename='C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//Int3480_2.png',
        color=inferno(100))

Rp1_m <- p[rownames(p) %in% Rp1,]
pheatmap(Rp1_m, 
         cluster_cols=FALSE, 
         show_rownames=FALSE, 
         show_colnames=FALSE,
         annotation_row=row_annotation, 
         annotation_colors=my_col, 
         cellheight=1.25, cellwidth=5, 
         legend=FALSE,
         annotation_legend=FALSE, 
         annotation_names_row=FALSE,
         border_color=NA,
         breaks=mat_breaks, 
         color=inferno(100),
         filename='C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//Rp1_2.png',
         treeheight_row=30
         )

RppM_m <- p[rownames(p) %in% RppM,]
pheatmap(RppM_m, 
         cluster_cols=FALSE, 
         show_rownames=FALSE, 
         show_colnames=FALSE,
         annotation_row=row_annotation, 
         annotation_colors=my_col, 
         cellheight=1.25, cellwidth=5, 
         legend=FALSE,
         annotation_legend=FALSE, 
         annotation_names_row=FALSE,
         annotation_names_col=FALSE,
         border_color=NA,
         annotation_col=col_annotation, #as the leftmost plot it gets the color annotation 
         breaks=mat_breaks, 
         color=inferno(100), 
        filename='C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//RppM_2.png',
         treeheight_row=30
         )


RppC_m <- p[rownames(p) %in% RppC,]
pheatmap(RppC_m, 
         cluster_cols=FALSE, 
         show_rownames=FALSE, 
         annotation_row=row_annotation, 
         annotation_colors=my_col, 
         cellheight=1.25, cellwidth=5, 
         legend=FALSE,
         annotation_legend=FALSE, 
         annotation_names_row=FALSE,
         breaks=mat_breaks, 
         color=inferno(100),
         show_colnames = FALSE,
         treeheight_row=30,
         filename='C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//RppC_2.png',
         border_color=NA)

#using this for the legend 
pheatmap(RppC_m, 
         cluster_cols=FALSE, 
         show_rownames=FALSE, 
         annotation_row=row_annotation, 
         annotation_colors=my_col, 
         cellheight=1.25, cellwidth=5, 
         legend=TRUE,
         annotation_legend=TRUE, 
         annotation_names_row=FALSE,
         breaks=mat_breaks, 
         color=inferno(100),
         show_colnames = TRUE,
         treeheight_row=30,
         filename='C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//legend.png',
                 border_color=NA
         )
```

## Supplemental Fig 3

Goal: create a standard violin plot between accessions in leaf tissue
(which tissue shows overall highest NLR expression?)

``` r
y_label <- expression('log'[2]*'(TPM)')

# group non-hvNLRs for ease of color plotting, add subpopulation information for grouping, then filter to just middle leaf tissue 
sfig3<- all_tpm_avg %>% 
  mutate(Clade_adj2=recode(Clade_adj, 'Int3480'='Int3480', 
                              'RppM-like'='RppM-like', 
                              'RppC-like'='RppC-like', 
                              'Rp1-like'='Rp1-like', 
                            .default='non-hvNLR')) %>%
  merge(subpopulations) %>% 
  arrange(subpopulation) %>%
  filter(tissue=='middle')

#add sample size of #NLRs to label 
sample_size <- sfig3 %>% 
  group_by(accession) %>% 
  summarize(n=n()) %>% 
  mutate(label=paste(accession, '\nn=', as.character(n), sep='')) %>% 
  ungroup() 

#factor and order
sfig3 <- merge(sfig3, sample_size)
subpop_order <- sfig3 %>% arrange(subpopulation) %>% pull(label) %>% unique()
sfig3$label <- factor(sfig3$label, levels=subpop_order)
sfig3$Clade_adj2 <- factor(sfig3$Clade_adj2, levels=c('non-hvNLR', 'Int3480', 'RppM-like', 'RppC-like', 'Rp1-like'))

#separate out tropical from the rest 
p1 <- sfig3 %>% filter(subpopulation != 'Tropical')
p2 <- sfig3 %>% filter(subpopulation == 'Tropical')
```

ggbeeswarm behaves poorly with faceting, so add labels manually

``` r
non_tropical <- ggplot(p1 %>% arrange(Clade_adj2), aes(x=label, y=log2_TPM))+ 
  geom_violin(trim=T,
              alpha = 0.4,
              scale='count')+
  geom_beeswarm(aes(color=Clade_adj2), 
                corral='random', 
                corral.width=0.8,
                cex=1,
                method='hex',
                priority='density',
                size=0.4)+
  scale_color_manual(values=c('Int3480'='#fee090', 
                              'RppM-like'='#d73027', 
                              'RppC-like'='#fc8d59', 
                              'Rp1-like'='#4575b4', 
                            'non-hvNLR'='grey'))+
  ylab(y_label)+
  xlab('')+
  theme(legend.title=element_blank(),
        legend.position='none',
        text = element_text(size=10))+
  ylim(0,10)+
  theme_classic() 

tropical <- ggplot(p2 %>% arrange(Clade_adj2), aes(x=label, y=log2_TPM))+ 
  geom_violin(trim=T,
              alpha = 0.4,
              scale='count')+
  geom_beeswarm(aes(color=Clade_adj2), 
                corral='random', 
                corral.width = 0.8,
                cex=1,
                method='hex',
                priority='density',
                size=0.4)+
  scale_color_manual(values=c('Int3480'='#fee090', 
                              'RppM-like'='#d73027', 
                              'RppC-like'='#fc8d59', 
                              'Rp1-like'='#4575b4', 
                            'non-hvNLR'='grey'))+
  ylab(y_label)+
  xlab('')+
  theme(legend.title=element_blank(),
        legend.position='none',
        text = element_text(size=10))+
  theme_classic() +
  ylim(0,10)

middle <- non_tropical + tropical + plot_layout(ncol=1, guides='collect')
ggsave(middle, filename='C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E16//figure panels//middle.png', width=250, height=150, units='mm')
```

    ## Warning: In `position_beeswarm`, method `hex` discretizes the data axis (a.k.a the
    ## continuous or non-grouped axis).
    ## This may result in changes to the position of the points along that axis,
    ## proportional to the value of `cex`.
    ## This warning is displayed once per session.

    ## Warning: Removed 524 rows containing missing values (`geom_point()`).

``` r
sfig3 %>% group_by(subpopulation, accession) %>% summarize(n=n())
```

    ## `summarise()` has grouped output by 'subpopulation'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 26 × 3
    ## # Groups:   subpopulation [6]
    ##    subpopulation   accession     n
    ##    <fct>           <chr>     <int>
    ##  1 Stiff stalk     B73         119
    ##  2 Non-stiff stalk B97         125
    ##  3 Non-stiff stalk KY21        126
    ##  4 Non-stiff stalk M162W       121
    ##  5 Non-stiff stalk MS71        122
    ##  6 Non-stiff stalk OH43        130
    ##  7 Non-stiff stalk OH7B        121
    ##  8 Mixed           M37W        134
    ##  9 Mixed           MO18W       117
    ## 10 Mixed           TX303       121
    ## # ℹ 16 more rows