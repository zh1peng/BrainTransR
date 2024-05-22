library(BrainEnrich)
library(dplyr)
brain_data=get_brainExample(type='PC1') %>% 
            filter(stringr::str_detect(Region, '^L_'))%>%
            tibble::column_to_rownames('Region')
gene_data=get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')

geneList.true=corr_brain_gene(gene_data, brain_data, method = 'pearson')  
geneSetList=get_geneSetList(type='CellType',parameter = 'Seidlitz2020')

geneSetList.filtered=filter_geneSetList(rownames(geneList.true), geneSetList, 20, 200)
gs_score.true=aggregate_geneSetList(geneList.true,geneSetList.filtered,  method = 'mean')

# resample gene
geneList.null=resample_gene(geneList.true, 1000)
gs_score.null=aggregate_geneSetList(geneList.null, geneSetList.filtered, method = 'mean')
pvals=caculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

# resample geneSetList matching coexpression
sampled_geneSetList=resample_geneSetList_matching_coexp(gene_data, geneSetList.filtered, tol = 0.1, max_iter = 1000000, n_perm = 5)
gs_score.null=aggregate_geneSetList_matching_coexp(geneList.true,geneSetList.filtered,sampled_geneSetList, method = 'mean')
pvals=caculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

# spin brain

