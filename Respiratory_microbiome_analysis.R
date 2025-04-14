### Anterior Nares and Nasopharynx Microbiome

library(phyloseq)
library(dplyr)
library(ggplot2)
library(microbiome)
library(microViz)


## phyloseq obj AN and NP samples

tax_mat <- taxonomy_RespSample %>%
  tibble::column_to_rownames("otu")

samples_df <- metadata_RespSample %>%
  tibble::column_to_rownames("SampleID")

otu_mat <- otu_RespSample %>%
  tibble::column_to_rownames("otu") 


otu_mat <- as.matrix(otu_mat)

tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

TAX = tax_table(tax_mat)

samples = sample_data(samples_df)

Resp_ALL <- phyloseq(OTU, TAX, samples)   ## 145 samples
Resp_ALL


### beta diversity comparison

library(MicrobiotaProcess)

Resp_ALL_MPSE <- Resp_ALL %>% as.MPSE()
Resp_ALL_MPSE

Resp_ALL_MPSE %<>% 
  mp_decostand(.abundance=Abundance)
Resp_ALL_MPSE

Resp_ALL_MPSE %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray") 
Resp_ALL_MPSE

Resp_ALL_MPSE %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")


Resp_ALL_MPSE %>% 
  mp_anosim(.abundance=hellinger, .group=sample_type, action="get") ## ANOSIM


Resp_ALL_MPSE %<>%
  mp_adonis(.abundance=hellinger, .formula=~sample_type, distmethod="bray", permutations=9999, action="add")
Resp_ALL_MPSE %>% mp_extract_internal_attr(name=adonis) ## ADONIS


beta <- Resp_ALL_MPSE %>%
  mp_plot_ord(
    .ord = pcoa, 
    .group = sample_type,
    .starshape = sample_type,
    .color = sample_type,
    .size = 2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE 
  ) +
  scale_fill_manual(values=c('orange2','lightblue')) +
  scale_color_manual(values=c('orange2','lightblue')) 
beta



## plot
library(ggpubr)

data_AlphaPlot_clean$higher = ifelse(data_AlphaPlot_clean$Anterior_Nares > data_AlphaPlot_clean$Nasopharynx, 'AN','NP')

alpha_comp <-ggpaired(data_AlphaPlot_clean,
                      cond1 = "Anterior_Nares", cond2 = "Nasopharynx",line.color = 'higher',
                      fill = "condition", palette = c('orange2','lightblue'), xlab = "", ylab = "",title = 'Alpha diversity (Shannon Index)')+ facet_wrap(~ Time)

alpha_comp

plotalpha_comp <- alpha_comp + theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),legend.title = element_text(size = 10),
                                     legend.text = element_text(size = 10),legend.position = 'right')
plotalpha_comp




