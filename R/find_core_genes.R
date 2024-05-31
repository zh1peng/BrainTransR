#' Find Core Genes Influencing Aggregated Score
#'
#' This function performs a Leave-One-Out (LOO) analysis on gene sets to determine core genes
#' that significantly influence the aggregated score. It utilizes parallel processing
#' to enhance computation efficiency.
#'
#' @param geneList A matrix of genes by models, each column representing a model.
#' @param geneSetList A list of gene sets, each containing names of genes.
#' @param method The aggregation method used to compute the scores.
#' @param n_cores The number of cores to use for parallel processing; defaults to 1.
#'                 Uses all available cores minus one if set to 0.
#' @param threshold_type The method to determine significance ('sd' for standard deviation, 'percentile' for percentile threshold).
#' @param threshold Numeric value specifying the threshold level; meaning depends on `threshold_type`.
#' @return A list of core genes for each gene set.
#' @export
find_core_genes <- function(geneList, geneSetList, method, n_cores = 1, threshold_type = c("sd", "percentile"), threshold = 1) {
  require(pbapply)
  require(parallel)

  # Validate the threshold_type argument
  threshold_type <- match.arg(threshold_type)


  if (threshold_type == "sd") {
    if (threshold < 0 || threshold > 3) {
      stop("For 'sd', threshold should be between 0 and 3.")
    }
  } else if (threshold_type == "percentile") {
    if (threshold < 1 || threshold > 99) {
      stop("For 'percentile', threshold should be a multiple of 10 and between 1 and 99.")
    }
  }

  # Determine the number of cores to use
  if (n_cores == 0 | n_cores > detectCores() - 1) {
    n_cores <- detectCores() - 1
  }

  # Initialize a cluster of workers
  cl <- makeCluster(n_cores)

  # Export necessary variables to the cluster
  clusterExport(cl, varlist = c("geneList", "aggregate_geneSet", "method"), envir = environment())

  # Parallelize the processing using pblapply for progress bar
  loo_changes <- pblapply(seq_along(geneSetList), function(i) {
    gs <- geneSetList[[i]]
    full_score <- aggregate_geneSet(geneList, gs, method)

    # Perform Leave-One-Out Analysis for each gene in the gene set
    loo_results <- sapply(gs, function(g) {
      modified_geneSet <- setdiff(gs, g)
      loo_score <- aggregate_geneSet(geneList, modified_geneSet, method)
      abs(full_score - loo_score)
    }, USE.NAMES = TRUE)

    names(loo_results) <- gs
    loo_results
  }, cl = cl)

  # Stop the cluster after processing
  stopCluster(cl)

  # Apply threshold to identify core genes
  if (threshold_type == "sd") {
    core_genes <- lapply(loo_changes, identify_core_genes_sd, k = threshold)
  } else if (threshold_type == "percentile") {
    core_genes <- lapply(loo_changes, identify_core_genes_percentile, percentile = threshold)
  }
  names(core_genes) <- names(geneSetList)
  return(core_genes)
}

#' Identify core genes based on a percentile threshold
#'
#' @param changes Named vector of changes from LOO analysis.
#' @param percentile Numeric value indicating the percentile to use as the threshold.
#' @return Vector of core genes.
identify_core_genes_percentile <- function(changes, percentile = 90) {
  threshold <- quantile(changes, probs = percentile / 100)
  core_genes <- names(changes[changes > threshold])
  return(core_genes)
}

#' Identify core genes based on a standard deviation threshold
#'
#' @param changes Named vector of changes from LOO analysis.
#' @param k Numeric value indicating how many standard deviations to use as the threshold.
#' @return Vector of core genes.
identify_core_genes_sd <- function(changes, k = 1) {
  mean_change <- mean(changes)
  sd_change <- sd(changes)
  threshold <- mean_change + k * sd_change
  core_genes <- names(changes[changes > threshold])
  return(core_genes)
}
