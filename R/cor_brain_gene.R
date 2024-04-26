
#' Calculate correlations or associations between gene and brain data
#'
#' @param gene_data A data frame or matrix of gene expression data.
#' @param brain_data A data frame or matrix of brain data.
#' @param method The method to be used for correlation/association. Can be 'pearson', 'spearman', 
#'               'pls1c', 'pls1w', or a custom function provided by the user.
#' @param r2z Logical, indicating whether to convert correlation coefficients to Fisher's Z scores.
#'            Only applicable to 'pearson'.
#'
#' @return A matrix with correlation or association coefficients between gene data and brain data.
#' @export
#' @examples
#' # Custom correlation function example
#' custom_cor <- function(gene_data, brain_data) {
#'   apply(gene_data, 2, function(g) apply(brain_data, 2, function(b) cor(g, b, method = "pearson")))
#' }
#' result <- corr_brain_gene(gene_data, brain_data, method = custom_cor)
corr_brain_gene <- function(gene_data, 
                            brain_data, 
                            method = c('pearson', 
                                       'spearman', 
                                       'pls1c', 
                                       'pls1w',
                                       'custom'), 
                            r2z = TRUE) {
  # Validate method
  if (is.function(method)) { # if method is a custom function
    assoc_func <- method
    method <- "custom"
  } else {
    assoc_func <- match.arg(method)
  }
  
  # Ensure rownames are matched for correlation calculation
  if (!identical(rownames(brain_data), rownames(gene_data))) {
    stop("Row names of brain_data and gene_data must match.")
  }
  
  # Initialize geneList
  geneList <- NULL
  
  # Define correlations or associations
  switch(assoc_func,
         pearson = {
           geneList <- cor(gene_data, brain_data, method = "pearson")
           if (r2z) geneList <- atanh(geneList)
         },
         spearman = {
           geneList <- cor(gene_data, brain_data, method = "spearman")
         },
         pls1c = {
           X <- apply(gene_data, 2, scale)
           y <- apply(brain_data, 2, scale)
           # do pls1 for each column
           pls1 <- apply(y,2,function(y_col) pls::plsr(y_col ~ X, 1))
           # extract the coefficients for the first component
           coef_list <- lapply(pls1,function(pls_i) drop(pls_i$coefficients))
           geneList=do.call(cbind,coef_list)
         },
         pls1w = {
            X <- apply(gene_data, 2, scale)
            y <- apply(brain_data, 2, scale)
            # do pls1 for each column
            pls1 <- apply(y,2,function(y_col) plsr(y_col ~ X, 1))
            # extract the coefficients for the first component
            weights_list <- lapply(pls1,function(pls_i) drop(loading.weights(pls_i)))
            geneList=do.call(cbind,weights_list)
         },
         custom = {
           geneList <- assoc_func(gene_data, brain_data)
         },
         {
           stop("Unsupported correlation method provided.")
         }
  )
  
  # Add attributes to geneList
  attr(geneList, "is_fisherz") <- if (assoc_func %in% c("pearson") && r2z) TRUE else NULL
  attr(geneList, "n.region") <- nrow(gene_data)
  attr(geneList, "cor_type") <- assoc_func
  
  return(geneList)
}



  




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

