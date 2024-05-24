
# This file contains miscellaneous functions that are used in the main scripts
caculate_pvals <- function(statList.true, statList.null,method=c('standard','split_pos_neg')) {
  
  method=match.arg(method)
  if (!identical(names(statList.true),names(statList.null))){
    stop('statList.true and statList.null are not matched')
  }

  if (method=='standard'){
  pvals=purrr::map2(statList.true, statList.null, function(true_stat, null_stat) {
    sum(abs(null_stat) >= abs(true_stat)) / (length(null_stat) + 1)
  })} else if (method=='split_pos_neg'){
  pvals=purrr::map2(statList.true, statList.null, function(true_stat, null_stat) {
    ifelse(true_stat >= 0,
           (sum(null_stat >= true_stat) + 1) / (sum(null_stat >= 0) + 1),
           (sum(null_stat <= true_stat) + 1) / (sum(null_stat < 0) + 1))
  })  
  }
  return(pvals)
}

calculate_qvalue <- function(pvals) {
    if (length(pvals) == 0)
        return(numeric(0))

    qobj <- tryCatch(qvalue::qvalue(pvals, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
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
    user_input <- readline(prompt = sprintf("%. Do you want to continue? (Y/N): ", msg))
    
    if (toupper(user_input) == "Y") {
      return(TRUE)
    } else if (toupper(user_input) == "N") {
      return(FALSE)
    } else {
      cat("Invalid input. Please enter 'Y' for Yes or 'N' for No.\n")
    }
  }
}