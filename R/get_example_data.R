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
load_example_data <- function(type=c('PC1','Cushing')) {
  file_path <- system.file("extdata", filename, package = "yourPackageName")
  
  if (file_path == "") {
    stop(sprintf('File %s does not exist in the package', filename))
  }
  
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  return(data)
}