load_permid <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainInfo/perm_id_weights',
                             atlas=c('desikan',
                                    'schaefer100',
                                    'schaefer200',
                                    'schaefer300'),
                             type=c('spin_brain','random_brain'),
                             perm.n=5000){
# load pre-generated perm.id
atlas=match.arg(atlas)
type=match.arg(type)

if (perm.n>10000){stop('>10k perm.id is supported')} 

 perm.id.file=sprintf('%s/%s_%s_perm_id.rds',data_path, atlas,type)
 
 # check if perm.id.file exist
  if (!file.exists(perm.id.file)){
    stop(sprintf('perm.id file %s does not exist',perm.id.file))
  }

 perm.id=readRDS(perm.id.file)   
  return(perm.id[,1:perm.n])
}






#' Generate Null Brain Data
#'
#' This function generates null brain datasets based on permutations provided in perm.id.
#' It rearranges the brain data according to the permutations and outputs the shuffled datasets.
#'
#' @param brain_data A matrix representing brain data, where each row corresponds to a region.
#' @param perm_id A matrix of permutations, where each column represents a permutation and each row corresponds to an index in `brain_data`.
#' 
#' @return A matrix of null brain data with the same dimensions as `brain_data` but with permuted rows according to `perm_id`.
#' @export
#'
#' @examples
#' brain_data <- matrix(rnorm(100), nrow=10)
#' perm_id <- matrix(sample(1:10), nrow=10, ncol=5)
#' null_data <- generate_null_brain_data(brain_data, perm_id)
generate_null_brain_data <- function(brain_data, perm_id) {
  # Check for duplicates in perm_id
  if (any(duplicated(perm_id))) {
    stop("Duplicate entries found in perm_id.")
  }

  # Define dimensions
  region_n <- nrow(perm_id)
  perm_n <- ncol(perm_id)

  # Check if dimensions match
  if (nrow(brain_data) != region_n) {
    stop("The number of regions in brain_data and perm_id do not match.")
  }

  # Generate null brain data
  null_brain_data <- sapply(1:perm_n, function(idx) {
    brain_data[perm_id[, idx], , drop = FALSE]
  })

  # Set row and column names
  rownames(null_brain_data) <- rownames(brain_data)
  colnames(null_brain_data) <- paste0("null_", 1:perm_n)

  return(null_brain_data)
}
