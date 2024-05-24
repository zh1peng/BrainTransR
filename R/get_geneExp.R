#' Get Gene Expression Data
#'
#' This function retrieves gene expression data based on specified parameters.
#' The data is obtained from the ENIGMA-TOOLBOX. Please cite the ENIGMA-TOOLBOX:
#' Larivière, S., Paquola, C., Park, B. Y., Royer, J., Wang, Y., Benkarim, O., ... & Bernhardt, B. C. (2021).
#' The ENIGMA Toolbox: multiscale neural contextualization of multisite neuroimaging datasets. Nature Methods, 18(7), 698-700.
#'
#' @param atlas Character. The brain atlas to use. One of 'desikan', 'schaefer100', or 'schaefer200'.
#' @param rdonor Character. The donor region to use. One of 'r0.2', 'r0.4', or 'r0.6'.
#' @param hem Character. The hemisphere to use. One of 'L', 'R', or 'B'.
#'
#' @return A matrix containing the gene expression data.
#' @export
#' @importFrom dplyr filter
#' @importFrom tibble column_to_rownames
#'
#' @examples
#' gene_exp_matrix <- get_geneExp(atlas = 'desikan', rdonor = 'r0.4', hem = 'L')
get_geneExp <- function(atlas = c('desikan', 'schaefer100', 'schaefer200', 'schaefer300'),
                        rdonor = c('r0.2', 'r0.4', 'r0.6'),
                        hem = c('L','R','B')) {
  
  atlas <- match.arg(atlas)
  rdonor <- match.arg(rdonor)
  hem <- match.arg(hem)

  if (atlas == 'desikan') {
    search_pattern <- sprintf('^%s_', hem)
  } else if (atlas %in% c('schaefer100', 'schaefer200', 'schaefer300')) {
    search_pattern <- sprintf('%sH', hem)
  }

  # Define the path to the gene expression data
  data_path <- system.file("extdata", package = "BrainEnrich")
  GeneExpCSV <- file.path(data_path, sprintf('geneExp/%s_%s.csv', atlas, rdonor))
  
  # Check if the CSV file exists, otherwise stop
  if (!file.exists(GeneExpCSV)) {
    stop(sprintf('GeneExp file %s does not exist', GeneExpCSV))
  }

  # Read the CSV file
  gene.df <- read.csv(GeneExpCSV, stringsAsFactors = FALSE)
  
  # Filter based on hemisphere
  if (hem %in% c('L', 'R')) {
    gene.df <- filter(gene.df, grepl(search_pattern, Region))
  }
  
  # Filter complete cases and convert to matrix
  gene.df <- gene.df %>%
    filter(complete.cases(.)) %>%
    column_to_rownames('Region')
  
  gene.mx <- as.matrix(gene.df)

  return(gene.mx)
}
