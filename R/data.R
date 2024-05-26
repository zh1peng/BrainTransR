#' Annotation Data for SynGO
#'
#' This dataset contains annotation data fetched from the BrainEnrich package for SynGO.
#'
#' @format A data frame with various annotation data.
#' @source BrainEnrich::get_annoData(type='SynGO')
"anno_data_synGO"

#' Brain Data Desikan Left Hemisphere PC1
#'
#' This dataset contains PC1 data filtered for regions starting with 'L_' in the Desikan atlas.
#'
#' @format A data frame with rows as regions and columns as PC1 data.
#' @source read.csv('data-raw/desikan_PC1_data.csv')
"brain_data_dk_lh_PC1"

#' Brain Data Desikan Left Hemisphere Cushing
#'
#' This dataset contains Cushing's data, filtered and reordered based on PC1 data for regions starting with 'L_' in the Desikan atlas.
#'
#' @format A data frame with rows as regions and columns as Cohen's d values.
#' @source read.csv('data-raw/desikan_cushing_data.csv')
"brain_data_dk_lh_Cushing"

#' Gene Expression Data Desikan Left Hemisphere
#'
#' This dataset contains gene expression data fetched from the BrainEnrich package for the Desikan atlas, left hemisphere.
#'
#' @format A data frame with various gene expression data.
#' @source BrainEnrich::get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')
"gene_data_dk_lh"

#' Desikan Centroid Coordinates for Left Hemisphere
#'
#' This dataset contains the centroid coordinates for the Desikan atlas regions in the left hemisphere.
#'
#' @format A data frame with rows as regions and columns as coordinates (x, y, z).
#' @source read.csv('data-raw/desikan_centroid.csv')
"coord_dk_lh"

#' Permutation IDs for Desikan Left Hemisphere
#'
#' This dataset contains permutation IDs generated for the Desikan spin brain for the left hemisphere.
#'
#' @format A matrix with permutation IDs for regions.
#' @source rotate_parcellation(coord.l = coord_dk_lh, nrot = 1000, seed = 2024)
"perm_id_dk_lh"
