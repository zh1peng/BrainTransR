#' Resample Gene List
#'
#' Generates a set of permuted gene lists from the original gene list, 
#' ensuring uniqueness in the permuted sets.
#'
#' @param geneList.true A matrix of gene expression data.
#' @param n_perm Number of permutations to generate.
#' @return A matrix of permuted gene lists.
#' @examples
#' gene_data <- matrix(rnorm(100), nrow=10)
#' resampled_gene_list <- resampling_geneList(gene_data, n_perm=100)
#' @export
resample_geneList <- function(geneList.true, n_perm = 5000) {

  # Check if the number of permutations requested exceeds 10,000. If it does, stop the function and display an error message.
  if (n_perm > 10000) {
    stop("Up to 10k permutations supported")
  }
  
  # Create a matrix 'geneList.null' by replicating gene lists. This is done 'n_perm + 100' times. Each replication involves shuffling the original gene list without replacement.
  # 'simplify = TRUE' ensures the result is a matrix.
  geneList.null <- replicate(n = n_perm + 100, expr = {
    sample(geneList.true, size = nrow(geneList.true), replace = FALSE)
  }, simplify = TRUE)
  
  # Remove duplicate columns from the transposed 'geneList.null' matrix. This is to ensure that all permutations are unique.
  geneList.null <- geneList.null[, !duplicated(t(geneList.null))]
  
  # Trim the matrix to keep only the first 'n_perm' columns, in case more than 'n_perm' unique permutations were generated initially.
  geneList.null <- geneList.null[, 1:n_perm]
  
  # Set the row names of 'geneList.null' to match those of 'geneList.true' to maintain consistency in identifiers.
  rownames(geneList.null) <- rownames(geneList.true)
  
  # Create column names for 'geneList.null' that indicate these are null permutations, e.g., 'null_1', 'null_2', ..., 'null_n_perm'.
  colnames(geneList.null) <- paste0("null_", seq_len(n_perm))

  # Copy attributes from 'geneList.true' to 'geneList.null' to retain metadata, such as dimension names, class, or other user-defined attributes.
  attr(geneList.null,'is_fisherz') <- attr(geneList.true,'is_fisherz')
  attr(geneList.null,'n.region') <- attr(geneList.true,'n.region')

  # Return the matrix of null permutations.
  return(geneList.null)
}


#' Resample Gene Sets with Specified Constraints
#'
#' This function resamples gene sets based on specific constraints like matching
#' co-expression patterns. The methodology implemented is informed by Wei et al. (2022) on statistical testing in transcriptomic-neuroimaging studies.
#' It is important to note that restricting null models to a subset of genes can be problematic. The empirical statistics sampled from 
#' the full gene pool differ from those derived from a restricted pool. Therefore, usage of this approach should be 
#' with caution.
#'
#' @param gene_data A matrix or data frame representing gene expression data.
#' @param geneSetList A list of gene sets to be resampled.
#' @param tol A numeric value indicating the tolerance for matching co-expression patterns (default = 0.01).
#' @param max_iter An integer indicating the maximum number of iterations for the sampling process (default = 1000000).
#' @param n_perm An integer indicating the number of permutations to generate (default = 5000).
#' @return A list of resampled gene sets based on the specified constraints.
#' @examples
#' # Assuming gene_data and geneSetList are predefined:
#' resampled_sets_coexp <- resample_geneSetList_coexp_matched(gene_data,
#'                                      geneSetList = geneSetList, tol = 0.01,
#'                                      max_iter = 1000000, n_perm = 5000)
#'
#' @references
#' Wei, Y., de Lange, S. C., Pijnenburg, R., Scholtens, L. H., Ardesch, D. J., Watanabe, K., Posthuma, D., & van den Heuvel, M. P. (2022).
#' Statistical testing in transcriptomic-neuroimaging studies: A how-to and evaluation of methods assessing spatial and gene specificity.
#' Human Brain Mapping, 43(3), 885–901. \url{https://doi.org/10.1002/hbm.25711}
resample_geneSetList_coexp_matched <- function(gene_data, geneSetList, tol = 0.01, max_iter = 1000000, n_perm = 5000) {
  # Calculate the co-expression matrix
  coexp_matrix <- cor(gene_data)
  
  # Sample each gene set to match co-expression patterns
  sampled_geneSetList <- lapply(geneSetList, function(gs) {
    message(sprintf('Sampling gene set %s\n', gs))
    sample_gs_matching_coexp(gs = gs, coexp_matrix = coexp_matrix, tol = tol, max_iter = max_iter, n_target = n_perm)
  })
  
  return(sampled_geneSetList)
}



#' Sample Gene Sets Matching Co-Expression
#'
#' This function samples gene sets that closely match the co-expression profile of a target gene set.
#' The methodology implemented is informed by Wei et al. (2022) on statistical testing in transcriptomic-neuroimaging studies.
#'
#' @param gs The target gene set for which similar co-expression profiles are sought.
#' @param coexp_matrix A co-expression matrix, typically calculated as the correlation matrix of gene expression data.
#' @param tol Tolerance for the difference between the co-expression of the target and the sampled gene sets.
#' @param max_iter Maximum number of iterations to attempt finding matches.
#' @param n_target Number of gene sets to sample.
#' @return A list of gene sets that closely match the target gene set's co-expression profile.
#' @references
#' Wei, Y., de Lange, S. C., Pijnenburg, R., Scholtens, L. H., Ardesch, D. J., Watanabe, K., Posthuma, D., & van den Heuvel, M. P. (2022).
#' Statistical testing in transcriptomic-neuroimaging studies: A how-to and evaluation of methods assessing spatial and gene specificity.
#' Human Brain Mapping, 43(3), 885–901. \url{https://doi.org/10.1002/hbm.25711}

#' Assuming 'coexp_matrix' is predefined:
#' target_gene_set <- c("Gene1", "Gene2", "Gene3")
#' similar_gene_sets <- sample_gs_matching_coexp(target_gene_set, coexp_matrix)
sample_gs_matching_coexp <- function(gs, coexp_matrix, tol = 0.01, max_iter = 1000000, n_target = 5000) {
  # Calculate the mean co-expression value for the lower triangle of the co-expression matrix of the target gene set
  gs_coexp_matrix <- coexp_matrix[gs, gs]
  gs_coexp_lower <- gs_coexp_matrix[lower.tri(gs_coexp_matrix)]
  target_coexp <- mean(gs_coexp_lower)

  # Initialize a list to hold sampled gene sets that meet the criteria
  sampled_gs <- list()

  for (i in 1:max_iter) {
    print(sprintf('Iteration %d', i))
    sampled_genes <- sample(colnames(coexp_matrix), length(gs), replace = FALSE)
    sampled_coexp_matrix <- coexp_matrix[sampled_genes, sampled_genes]
    sampled_coexp_lower <- sampled_coexp_matrix[lower.tri(sampled_coexp_matrix)]
    sampled_coexp_avg <- mean(sampled_coexp_lower)

    if (abs(sampled_coexp_avg - target_coexp) < tol) {
      sampled_gs[[length(sampled_gs) + 1]] <- sampled_genes
    }

    # Remove duplicates and check if the target number of sets has been reached
    sampled_gs <- unique(sampled_gs)
    if (length(sampled_gs) == n_target) {
      break
    }
  }

  if (length(sampled_gs) < n_target) {
    stop('n_target not reached, consider increasing max_iter or decreasing tol')
  }

  return(sampled_gs)
}



