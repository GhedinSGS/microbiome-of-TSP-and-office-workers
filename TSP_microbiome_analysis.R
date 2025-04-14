### TSP microbiome

library(phyloseq)
library(dplyr)
library(ggplot2)
library(microbiome)
library(microViz)

## create phyloseq obj.

tax_mat <- taxonomy_TSP %>%
  tibble::column_to_rownames("otu")

samples_df <- metadata_TSP %>%
  tibble::column_to_rownames("SampleID")

otu_mat <- otu_table_TSP %>%
  tibble::column_to_rownames("otu") 


otu_mat <- as.matrix(otu_mat)

tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

TAX = tax_table(tax_mat)

samples = sample_data(samples_df)

TSP_ALL <- phyloseq(OTU, TAX, samples)   ## all TSP samples
TSP_ALL






## select INDOOR samples
TSP_Indoor <- ps_filter(TSP_ALL, TSP_type == 'Indoor')
TSP_Indoor

## create a MPSE file
library(MicrobiotaProcess)

TSP_Indoor <- TSP_Indoor %>% as.MPSE()
TSP_Indoor

TSP_Indoor %<>% 
  mp_cal_alpha(.abundance=Abundance) ## with rarefraction
TSP_Indoor

alpha_Indoor <- TSP_Indoor %>% 
  mp_plot_alpha(
    .group=Groups, 
    .alpha=c('Shannon'),
    test = "wilcox.test",
    size = 0.2,
    text = 4)

alpha_Indoor





### select OUTDOOR samples
TSP_Outdoor <- ps_filter(TSP_ALL, TSP_type == 'Outdoor')
TSP_Outdoor

## alpha diversity

TSP_Outdoor <- TSP_Outdoor %>% as.MPSE()
TSP_Outdoor

TSP_Outdoor %<>% 
  mp_cal_alpha(.abundance=Abundance) ## with rarefraction
TSP_Outdoor

alpha_Outdoor <- TSP_Outdoor %>% 
  mp_plot_alpha(
    .group=Groups, 
    .alpha=c('Shannon'),
    test = "wilcox.test",
    size = 0.2,
    text = 4) 

alpha_Outdoor






### OUTDOOR vs INDOOR

TSP_IndoorOutdoor <- TSP_ALL %>% as.MPSE()
TSP_IndoorOutdoor

TSP_IndoorOutdoor %<>% 
  mp_cal_alpha(.abundance=Abundance) ## with rarefraction
TSP_IndoorOutdoor

alpha_IndoorOutdoor <- TSP_IndoorOutdoor %>% 
  mp_plot_alpha(
    .group=TSP_type, 
    .alpha=c('Shannon'),
    test = "wilcox.test",
    size = 0.2,
    text = 4)  +
  scale_fill_manual(values=c('palegreen3', "lightgrey")) +
  scale_color_manual(values=c("palegreen3", "lightgrey"))

alpha_IndoorOutdoor


TSP_IndoorOutdoor %<>% 
  mp_decostand(.abundance=Abundance)
TSP_IndoorOutdoor

TSP_IndoorOutdoor %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
TSP_IndoorOutdoor

TSP_IndoorOutdoor %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")


TSP_IndoorOutdoor %>% 
  mp_anosim(.abundance=hellinger, .group=TSP_type, action="get") ## ANOSIM


TSP_IndoorOutdoor %<>%
  mp_adonis(.abundance=hellinger, .formula=~TSP_type, distmethod="bray", permutations=9999, action="add")
TSP_IndoorOutdoor %>% mp_extract_internal_attr(name=adonis) ## ADONIS R2


beta_IndoorOutdoor <- TSP_IndoorOutdoor %>%
  mp_plot_ord(
    .ord = pcoa, 
    .group = TSP_type,
    .starshape = TSP_type,
    .color = TSP_type,
    .size = 2.5,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse
  ) +
  scale_fill_manual(values=c("lightgrey",'palegreen3')) +
  scale_color_manual(values=c('lightgrey','palegreen3')) 

beta_IndoorOutdoor



### ALDEX2 analysis (Indoor vs Outdoor)
library(ALDEx2)

## select top 200 for the analysis
IndoorOutdoor_200 = top_taxa(TSP_ALL, n=200)

IndoorOutdoor_top = prune_taxa(IndoorOutdoor_200,TSP_ALL)
IndoorOutdoor_top

TSP_200 = as.data.frame(IndoorOutdoor_top@otu_table)

TSP_200 = tibble::rownames_to_column(TSP_200,'otu')
TSP_200 = merge(taxonomy_TSP,TSP_200,by = 'otu')
TSP_200 = TSP_200[,-1]
TSP_200 = tibble::column_to_rownames(TSP_200, 'Species')

TSP_200_condition = as.data.frame(t(TSP_200))
TSP_200_condition = tibble::rownames_to_column(TSP_200_condition,'SampleID')
TSP_200_condition = TSP_200_condition[,1:2]
TSP_200_condition = merge(TSP_200_condition,metadata_TSP,by = 'SampleID')
TSP_condition = TSP_200_condition[,4]

## analysis

prova.aldex <- aldex(TSP_200, TSP_condition, test="t", effect=TRUE, 
                     include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)

TSP_aldexResults = tibble::rownames_to_column(prova.aldex, 'Species')

## plot the results 

aldex.plot(prova.aldex) 

x.u <- aldex.clr(TSP_200, TSP_condition, verbose=FALSE)
par(col = "black")
plot <- aldex.plotFeature(x.u, 'Corynebacterium diphtheriae' ,pooledOnly=T,densityOnly = T)




