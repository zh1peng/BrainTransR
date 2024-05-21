#' Load example Data
#'
#' This function loads example data included with the package.
#'
#' @param filename Character. The name of the file to load.
#'
#' @return A data frame containing the example data.
#' @export
#' @importFrom utils read.csv
#'
#' @examples
#' example_data <- load_example_data('example.csv')
get_brainExample <- function(type=c('PC1','Cushing')) {
    
    type <- match.arg(type)
data_path <- system.file("extdata", package = "BrainEnrich")

filename <- switch(type,
                   PC1 = 'desikan_PC1_data.csv',
                   Cushing = 'desikan_cushing_data.csv',
                   stop("Invalid type")
                  )

  file_path <- file.path(data_path, sprintf('example/%s', filename))
  if (file_path == "") {
    stop(sprintf('File %s does not exist in the package', filename))
  }
  
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  return(data)
}