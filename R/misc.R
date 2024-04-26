

#' Calculate correlations or associations between gene and brain data
#'
#' @param statList.true A data frame or matrix of gene expression data.
#' @param statList.null A data frame or matrix of brain data.
#' @param method The method to be used for correlation/association. Can be 'pearson', 'spearman', 
#'               'pls1c', 'pls1w', or a custom function provided by the user.

#'
#' @return A matrix with correlation or association coefficients between gene data and brain data.
#' @export
#' @examples
#' # Custom correlation function example
#' custom_cor <- function(gene_data, brain_data) {
#'   apply(gene_data, 2, function(g) apply(brain_data, 2, function(b) cor(g, b, method = "pearson")))
#' }
#' result <- corr_brain_gene(gene_data, brain_data, method = custom_cor)
caculate_pvals <- function(statList.true, statList.null,method=c('standard','split_pos_neg')) {
  
  method=match.arg(method)
  if (!identical(names(statList.true),names(statList.null))){
    stop('statList.true and statList.null are not matched')
  }

  if (method=='standard'){
  pvals=purrr::map2(statList.true, statList.null, function(true_stat, null_stat) {
    sum(abs(null_stat) >= abs(true_stat)) / (length(null_stat) + 1)
  })} else if (method=='split_pos_neg'){
  pvals=purrr::map2(statList.true, statList.null, function(true_stat, null_stat) {
    ifelse(true_stat >= 0,
           (sum(null_stat >= true_stat) + 1) / (sum(null_stat >= 0) + 1),
           (sum(null_stat <= true_stat) + 1) / (sum(null_stat < 0) + 1))
  })  
  }
  return(pvals)
}




cor2p <- function(r,n){
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  return(p)
}
