library(BrainEnrich)
library(dplyr)
brain_data=get_brainExample(type='PC1') %>% 
            filter(stringr::str_detect(Region, '^L_'))%>%
            tibble::column_to_rownames('Region')
gene_data=get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')

geneList.true=corr_brain_gene(gene_data, brain_data, method = 'pearson')  
geneSetList=get_geneSetList(type='GO',parameter = 'CC')

geneSetList.filtered=filter_geneSetList(rownames(geneList.true), geneSetList, 20, 200)
gs_score.true=aggregate_geneSetList(geneSetList.filtered, geneList.true, method = 'mean')

# resampling gene
geneList.null=resample_gene(geneList.true, 1000)
gs_score.null=aggregate_geneSetList(geneSetList.filtered, geneList.null, method = 'mean')
pvals=caculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

# resampling gene with coexpression matched
geneList.null=resample_gene_coexp_matched(gene_data, geneSetList.filtered, tol = 0.01, max_iter = 1000000, n_perm = 1000)



