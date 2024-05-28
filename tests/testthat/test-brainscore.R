library(BrainEnrich)
data('brain_data_random')
data('gene_data_sample')
data('annoData_synGO')

gene_data <- gene_data_sample
brain_data <- brain_data_random
cor_method <- 'pearson'
aggre_method <- 'mean'
prefix <- 'SynGO_'
annoData <- annoData_synGO
res=brainscore(brain_data=brain_data, gene_data=gene_data, annoData,cor_method, aggre_method, prefix)




geneList <- corr_brain_gene(gene_data, brain_data, method = cor_method) 
geneSetList <- BrainEnrich:::get_geneSetList(annoData)
selected.gs <- filter_geneSetList(rownames(geneList), geneSetList, minGSSize = 10, maxGSSize = 200) 
gs.score <- aggregate_geneSetList(geneList, selected.gs, method = aggre_method, prefix = prefix)
res <- data.frame(gs.score)
rownames(res) <- colnames(brain_data)
colnames(res) <- names(selected.gs)

print(res)








