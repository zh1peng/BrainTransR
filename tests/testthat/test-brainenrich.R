library(BrainEnrich)
library(dplyr)
brain_data=get_brainExample(type='PC1') %>% 
            filter(stringr::str_detect(Region, '^L_'))%>%
            tibble::column_to_rownames('Region')
gene_data=get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')

geneList.true=corr_brain_gene(gene_data, brain_data, method = 'pearson')  
annoData=get_annoData(type='CellType',parameter = 'Seidlitz2020')
geneSetList=get_geneSetList(annoData)


selected.gs=filter_geneSetList(rownames(geneList.true), geneSetList, 20, 200)
gs_score.true=aggregate_geneSetList(geneList.true,selected.gs, n_cores = 0, prefix=NULL,  method = 'mean')

# resample gene
geneList.null=resample_gene(geneList.true, 1000)
gs_score.null=aggregate_geneSetList(geneList.null, selected.gs, method = 'mean')
pvals=calculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

# resample geneSetList matching coexpression
sampled_gs=resample_geneSetList_matching_coexp(gene_data, selected.gs, tol = 0.1, max_iter = 1000000, n_perm = 100)
gs_score.null=aggregate_geneSetList_matching_coexp(geneList.true,selected.gs,sampled_gs, method = 'mean')
pvals=calculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

# spin brain provide coord.l
perm_id=rotate_parcellation(coord.l = , nrot = 1000)
null_brain_data=generate_null_brain_data(brain_data, perm_id)
geneList.null=corr_brain_gene(gene_data, null_brain_data, method = 'pearson')
gs_score.null=aggregate_geneSetList(geneList.null, selected.gs, method = 'mean')
pvals=caculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

pvals.adj <- p.adjust(pvals, method='fdr')
qvalues <- calculate_qvalue(unlist(pvals))
# prepare results for group level




gs.name=names(selected.gs)
Description <- names(selected.gs)


params <- list(pvalueCutoff = 0.05,
                nPerm = 1000,
                pAdjustMethod = 'BH',
                minGSSize = 20,
                maxGSSize = 200
                   )
 

res <- data.frame(
  ID = as.character(gs.name),
  Description = as.character(Description),
  setSize = sapply(selected.gs, length),
  gsScore= unlist(gs_score.true),
  pvalue = unlist(pvals),
  p.adjust=pvals.adj,
  qvalue = qvalues,
  stringsAsFactors = FALSE
)
rownames(res) <- NULL

attributes(geneList.true) <- NULL

res=new("gseaResult",
        result     = res,
        geneSets   = selected.gs,
        geneList   = geneList.true,
        params     = params,
        readable   = FALSE
        )
res@organism <- "UNKNOWN"
res@setType <- "UNKNOWN"
res@keytype <- "UNKNOWN"

library(enrichplot)
barplot(res)
cnetplot(res)
heatplot(res,showCategory = 3)
treeplot(res)

emapplot(res,showCategory = 3)