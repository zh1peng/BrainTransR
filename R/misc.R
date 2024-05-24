#' Calculate P-Values
#'
#' This function calculates p-values based on the provided true statistics and null statistics lists.
#' It supports two methods for p-value calculation: 'standard' and 'split_pos_neg'.
#'
#' @param statList.true A named list of true statistics.
#' @param statList.null A named list of null statistics corresponding to the true statistics.
#' @param method The method to be used for p-value calculation. Either 'standard' or 'split_pos_neg'. Default is 'standard'.
#' @return A list of calculated p-values.
#' @examples
#' # Example usage of calculate_pvals function
#' true_stats <- list(A = 2.5, B = -1.8, C = 0.5)
#' null_stats <- list(A = rnorm(1000), B = rnorm(1000), C = rnorm(1000))
#' pvals <- calculate_pvals(true_stats, null_stats, method = 'standard')
#' print(pvals)
#' @importFrom purrr map2
#' @export
calculate_pvals <- function(statList.true, statList.null, method = c('standard', 'split_pos_neg')) {
  method <- match.arg(method)
  
  if (!identical(names(statList.true), names(statList.null))) {
    stop('statList.true and statList.null are not matched')
  }
  
  pvals <- map2(statList.true, statList.null, function(true_stat, null_stat) {
    if (method == 'standard') {
      sum(abs(null_stat) >= abs(true_stat)) / (length(null_stat) + 1)
    } else if (method == 'split_pos_neg') {
      if (true_stat >= 0) {
        (sum(null_stat >= true_stat) + 1) / (sum(null_stat >= 0) + 1)
      } else {
        (sum(null_stat <= true_stat) + 1) / (sum(null_stat < 0) + 1)
      }
    }
  })
  
  return(pvals)
}


#' Calculate Q-Values
#'
#' This function calculates q-values from a vector of p-values using the `qvalue` package.
#' Q-values are used to control the false discovery rate (FDR) in multiple hypothesis testing.
#'
#' @param pvals A numeric vector of p-values.
#' @return A numeric vector of q-values corresponding to the input p-values, or NA if q-value calculation fails.
#' @examples
#' pvals <- runif(100)
#' qvals <- calculate_qvalue(pvals)
#' print(qvals)
#' @importFrom qvalue qvalue
#' @export
calculate_qvalue <- function(pvals) {
  if (length(pvals) == 0) {
    return(numeric(0))
  }
  
  qobj <- tryCatch(
    qvalue(pvals, lambda = 0.05, pi0.method = "bootstrap"),
    error = function(e) NULL
  )
  
  if (!is.null(qobj) && class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- rep(NA, length(pvals))
  }
  
  return(qvalues)
}


cor2p <- function(r,n){
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  return(p)
}



ask_user_continue <- function(msg) {
  repeat {
    user_input <- readline(prompt = sprintf("%s. Do you want to continue? (Y/N): ", msg))
    
    if (toupper(user_input) == "Y") {
      return(TRUE)
    } else if (toupper(user_input) == "N") {
      return(FALSE)
    } else {
      cat("Invalid input. Please enter 'Y' for Yes or 'N' for No.\n")
    }
  }
}