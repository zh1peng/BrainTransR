library(BrainEnrich)
library(dplyr)

data("gene_data_sample")
data("brain_data_PC1")
data("annoData_synGO")
data("coord_dk_lh")

gene_data=gene_data_sample
annoData=annoData_synGO
brain_data=brain_data_PC1
coord.l=coord_dk_lh



geneList.true=corr_brain_gene(gene_data, brain_data, method = 'pearson')  
geneSetList=get_geneSetList(annoData)

selected.gs=filter_geneSetList(rownames(geneList.true), geneSetList, 20, 200)
gs_score.true=aggregate_geneSetList(geneList.true,selected.gs, n_cores = 0, prefix=NULL,  method = 'mean')


core_genes=find_core_genes(geneList.true, selected.gs, method = 'mean', n_cores = 0, threshold_type = 'sd', threshold = 1)

# resample gene
geneList.null=resample_gene(geneList.true, 1000)
gs_score.null=aggregate_geneSetList(geneList.null, selected.gs, method = 'mean')
pvals=calculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

# resample geneSetList matching coexpression
sampled_gs=resample_geneSetList_matching_coexp(gene_data, selected.gs, tol = 0.1, max_iter = 1000000, n_perm = 100)
gs_score.null=aggregate_geneSetList_matching_coexp(geneList.true,selected.gs,sampled_gs, method = 'mean')
pvals=calculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

# spin brain provide coord.l

perm_id=rotate_parcellation(coord.l = coord.l, nrot = 5000, seed=2024)
null_brain_data=generate_null_brain_data(brain_data, perm_id)
geneList.null=corr_brain_gene(gene_data, null_brain_data, method = 'pearson')
gs_score.null=aggregate_geneSetList(geneList.null, selected.gs, method = 'mean')
pvals=calculate_pvals(gs_score.true, gs_score.null, method=c('standard'))

pvals.adj <- p.adjust(pvals, method='fdr')
qvalues <- calculate_qvalue(unlist(pvals))
# prepare results for group level


TERM2NAME=getFromNamespace("TERM2NAME", "DOSE")
calculate_qvalue=getFromNamespace("calculate_qvalue", "DOSE")

gs.name=names(selected.gs)
Description <-  TERM2NAME(gs.name, annoData)


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

res$core_enrichment <- sapply(core_genes, paste0, collapse='/')




# attr(geneList.true, "is_fisherz") <- NULL
# attr(geneList.true, "n.region") <- NULL
# attr(geneList.true, "cor_type") <- NULL

geneList.null=

as.data.frame(gs_score.null)
a=bind_rows(gs_score.null)

res.obj=new("gseaResult",
        result     = res,
        geneSets   = selected.gs,
        geneList   = geneList.true[,1],
        params     = params,
        readable   = TRUE
        )
res.obj@organism <- "Homo sapiens"
res.obj@setType <- "UNKNOWN"
res.obj@keytype <- "Symbol"


library(enrichplot)
res.obj@result%>%barplot(x='pvalue')

heatplot(res.obj,showCategory = 3)
treeplot(res.obj)
emapplot(res,showCategory = 3)

upsetplot(res.obj,n=5)

list2df=getFromNamespace('list2df','enrichplot')
extract_geneSets=getFromNamespace('extract_geneSets','enrichplot')
geneInCategory=getFromNamespace('geneInCategory','enrichplot')
x=res.obj
 geneSets <- extract_geneSets(x, 10)

cnetplot(res.obj)
barplot(res.obj)
dotplot(res)
upsetplot(res.obj)


library(DOSE)
data(geneList)
edo2 <- gseDO(geneList)

cnetplot(edo2)

de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)

edox <- setReadable(edo2, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, categorySize="pvalue")

extract_geneSets(edo2,2)
edo2.df=as.data.frame(edo2)

extract_geneSets(res.obj,2)