#' SynGO Annotation Data
#'
#' This dataset includes synaptic functionalities and gene ontology annotations from the SynGO project,
#' useful for studies focused on synaptic functions and disorders.
#'
#' @format A \code{\link[data.frame]{data.frame}}
#' @source SynGO project, sourced from \url{https://www.syngoportal.org/}
#' @examples
#' data(syngo_annoData)
#' summary(syngo_annoData)
"syngo_annoData"

#' Permutation ID Data for Desikan Spin Brain (Left Hemisphere)
#'
#' Permutation identifiers used for analyses involving spatial autocorrelation in neuroimaging data,
#' specifically tailored for the left hemisphere based on the Desikan atlas.
#'
#' @format An R object loaded via \code{\link[base]{readRDS}}
#' @source Data generated internally using Desikan atlas configurations.
#' @examples
#' data(perm_id)
#' length(perm_id)
"perm_id"

#' Principal Component 1 Data Filtered by Region (Left Hemisphere)
#'
#' The first principal component of brain data from the Desikan atlas,
#' filtered to include only regions starting with 'L_', indicative of the left hemisphere.
#'
#' @format A \code{\link[data.frame]{data.frame}} with regions as row names.
#' @source Processed based on PC1 analysis from MRI data of the Desikan atlas (left hemisphere).
#' @examples
#' data(PC1_data)
#' head(PC1_data)
"PC1_data"

#' Cushing's Disease Related Data (Left Hemisphere)
#'
#' Processed data from a study on Cushing's disease, focusing on Cohen's d values for brain regions,
#' specifically ordered and filtered for the left hemisphere using the Desikan atlas.
#'
#' @format A \code{\link[data.frame]{data.frame}} with regions as row names.
#' @source Adapted from clinical study data on Cushing's disease, aligned with Desikan atlas regions (left hemisphere).
#' @examples
#' data(Cushing_data)
#' plot(Cushing_data$cohend, type = 'h')
"Cushing_data"

#' Gene Expression Data from Desikan Atlas (Left Hemisphere)
#'
#' This dataset provides gene expression levels across various regions of the left hemisphere,
#' based on the Desikan atlas, aimed at neurogenetic studies.
#'
#' @format A \code{\link[data.frame]{data.frame}} with gene identifiers as row names.
#' @source Gene expression profiles derived from the Desikan atlas, release 'r0.6', specifically for the left hemisphere.
#' @examples
#' data(gene_data)
#' summary(gene_data$expression_level)
"gene_data"
