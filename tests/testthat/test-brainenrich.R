library(BrainEnrich)
library(dplyr)
brain_data=get_brainExample(type='PC1') %>% 
            filter(stringr::str_detect(Region, '^L_'))%>%
            tibble::column_to_rownames('Region')
gene_data=get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')

geneList.true=corr_brain_gene(gene_data, brain_data, method = 'pearson')  
geneSetList=get_geneSetList(type='GO',parameter = 'CC')

