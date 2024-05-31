#' SynGO Annotation Data
#'
#' This dataset contains annotation data fetched from the BrainEnrich package for SynGO.
#'
#' @format A data frame with various annotation data.
#' @source BrainEnrich::get_annoData(type = 'SynGO')
"annoData_synGO"

#' Brain Data PC1 for Left Hemisphere
#'
#' This dataset contains PC1 data filtered for regions starting with 'L_' in the Desikan atlas.
#'
#' @format A data frame with rows as regions and columns as PC1 data.
#' @source read.csv('data-raw/desikan_PC1_data.csv')
"brain_data_PC1"

#' Brain Data Cushing for Left Hemisphere
#'
#' This dataset contains Cushing's data, filtered and reordered based on PC1 data for regions starting with 'L_' in the Desikan atlas.
#'
#' @format A data frame with rows as regions and columns as Cohen's d values.
#' @source read.csv('data-raw/desikan_cushing_data.csv')
"brain_data_Cushing"

#' Gene Expression Data for Desikan Atlas Left Hemisphere
#'
#' This dataset contains gene expression data fetched from the BrainEnrich package for the Desikan atlas, left hemisphere.
#'
#' @format A data frame with various gene expression data.
#' @source BrainEnrich::get_geneExp(atlas = 'desikan', rdonor = 'r0.6', hem = 'L')
"gene_data_sample"

#' Desikan Centroid Coordinates for Left Hemisphere
#'
#' This dataset contains the centroid coordinates for the Desikan atlas regions in the left hemisphere.
#'
#' @format A data frame with rows as regions and columns as coordinates (x, y, z).
#' @source read.csv('data-raw/desikan_centroid.csv')
"coord_dk_lh"




#' Random Desikan region for testing brainscore function
#'
#' @format Simulated data frame for 34 regions by 100 subject.
#' @source read.csv('data-raw/desikan_centroid.csv')
"brain_data_random"
