#' Aggregate Gene List Based on a Gene Set
#'
#' This function aggregates a gene list based on a specified gene set using various methods.
#'
#' @param geneList A matrix representing named correlation coefficients, with each column representing a model.
#' @param geneSet A vector of gene names representing the gene set of interest.
#' @param method A character string specifying the aggregation method to be used, or a custom function provided by the user.

#' @return Returns a numeric vector of aggregated gene scores.
#' @examples
#' geneList <- matrix(rnorm(100), nrow = 10)
#' rownames(geneList) <- paste("Gene", 1:10, sep = "")
#' geneSet <- rownames(geneList)[1:5]
#' aggregate_geneSet(geneList, geneSet, method = "mean")
#' @export
#' 
#' 
aggregate_geneSet <- function(geneList,  # named correlation/coefficient(pls1) matrix
                              geneSet,   # one geneSet of interest
                              method = c('mean', 
                                         'median', 
                                         'meanabs', 
                                         'meansqr',
                                         'maxmean', 
                                         'sig_n',
                                         'sign_test',
                                         'rank_sum',
                                         'ks_orig',
                                         'ks_weighted',
                                         'ks_sum',
                                         'locfdr')) {

                              }