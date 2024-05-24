#' Calculate Brain Scores
#'
#' This function calculates brain scores based on the provided brain data frame, gene data, and gene set list.
#' It computes correlations between brain regions and genes, and aggregates these scores based on the specified methods.
#'
#' @param brain_df Data frame with subjects as rows and brain regions as columns.
#' @param gene_data Data frame with brain regions as rows and genes as columns.
#' @param geneSetList A list of gene sets to aggregate scores for.
#' @param cor_method Character string specifying the correlation method. One of 'pearson', 'spearman', 'pls1c', 'pls1w', or 'custom'.
#' @param aggre_method Character string specifying the aggregation method. One of "mean", "median", "meanabs", "meansqr",
#' "maxmean", "ks_orig", "ks_weighted", "ks_pos_neg_sum", "local_fdr", "sign_test", "rank_sum", or "custom".
#' @param prefix Optional character string to prefix the names of the aggregated scores.
#' @return A data frame with aggregated gene set scores.
#' @examples
#' \dontrun{
#' # Example usage of brainscore function
#' brain_df <- matrix(rnorm(100), nrow = 10)
#' gene_data <- matrix(rnorm(100), nrow = 10)
#' geneSetList <- list(set1 = c("gene1", "gene2"), set2 = c("gene3", "gene4"))
#' scores <- brainscore(brain_df, gene_data, geneSetList, cor_method = 'pearson', aggre_method = 'mean')
#' print(scores)
#' }
#' @importFrom stats cor
#' @export
brainscore <- function(brain_df, 
                       gene_data, 
                       geneSetList,
                       cor_method = c('pearson', 'spearman', 'pls1c', 'pls1w', 'custom'),
                       aggre_method = c("mean", "median", "meanabs", "meansqr", "maxmean", 
                                        "ks_orig", "ks_weighted", "ks_pos_neg_sum", 
                                        "local_fdr", "sign_test", "rank_sum", "custom"),
                       prefix = NULL) {

  cor_method <- match.arg(cor_method)
  aggre_method <- match.arg(aggre_method)

  brain_data <- t(brain_df)
  geneList <- corr_brain_gene(gene_data, brain_data, method = cor_method)  
  gs.score <- aggregate_geneSetList(geneList, geneSetList, method = aggre_method, prefix = prefix)
  df.score <- data.frame(gs.score)
  
  return(df.score)
}
