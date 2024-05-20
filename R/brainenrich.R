brainenrich <-function( gene_data, 
                        brain_data,
                        annoData,
                        cor_method=c('pearson', 'spearman', 'pls1c', 
                                    'pls1w', 'custom'),
                        aggre_method= c("mean", "median", "meanabs", "meansqr",
                                        "maxmean","ks_orig", "ks_weighted", "ks_pos_neg_sum", 
                                        "local_fdr", "sign_test","rank_sum", "custom"),
                        null_model=c('spin_brain','resample_gene','coexp_matched'),
                        n_perm=1000,
                        coord.l=NULL,
                        coord.r=NULL,
                        minGSSize=10,
                        maxGSSize=200){

geneList.true=corr_brain_gene(gene_data, brain_data, method = cor_method)                      
geneSetList=filter_geneSetList(names(geneList.true), geneSetList, minGSSize, maxGSSize)
gs_score.true=aggregate_geneSetList(geneSetList, geneList.true, method = aggre_method)

if (null_model == 'spin_brain'){
perm_id=rotate_parcellation(coord.l = coord.l, coord.r = coord.r, nrot = n_perm)
null_brain_data=generate_null_brain_data(brain_data, perm_id)
geneList.null=corr_brain_gene(gene_data, null_brain_data, method = cor_method)

} else if (null_model == 'resample_gene'){ 
geneList.null=resample_geneList(geneList.true, n_perm)

} else if (null_model == 'coexp_matched'){ 
geneList.null=resample_geneSetList_coexp_matched(gene_data, geneSetList, tol = 0.01, max_iter = 1000000, n_perm = n_perm)
}

gs_score.null=aggregate_geneSetList(geneSetList, geneList.null, method = aggre_method)

pvals=caculate_pvals(gs_score.true, gs_score.null, method=c('standard'))


}



                                  