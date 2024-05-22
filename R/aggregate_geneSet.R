#' Aggregate Gene Set Scores
#'
#' Function to aggregate geneList based on geneSet, evaluating one geneSet at a time.
#' The function supports multiple aggregation methods as specified by the user.
#'
#' @param geneList A matrix of genes by models, with each column representing a true or null model.
#' @param geneSet A vector containing names of genes in the gene set of interest.
#' @param method A character string specifying the method to use for aggregation. 
#'               Options include 'mean', 'median', 'meanabs', 'meansqr', 'maxmean', 
#'               'sig_n', 'sign_test', 'rank_sum', 'ks_orig', 'ks_weighted', 'ks_sum', 
#'               and 'locfdr'. Default is 'mean'.
#
#' @return Returns a numeric score based on the specified aggregation method.
#' @examples
#' geneList <- matrix(rnorm(100), ncol=10, dimnames=list(gene=paste("Gene", 1:10), model=1:10))
#' gene_set <- paste("Gene", sample(1:10, 5), sep="")
#' aggregate_geneSet(gene_data, gene_set, method = 'mean')
#' @export
#' 
#' 
aggregate_geneSet <- function(geneList, # named correlation/coefficient matrix
                              geneSet, # one geneSet of interest
                              method = c(
                                  "mean", "median", "meanabs", "meansqr",
                                  "maxmean",   "ks_orig", "ks_weighted", 
                                  "ks_pos_neg_sum", "local_fdr", 
                                  "sign_test", "rank_sum", "custom")) {
    if (is.function(method)) { # if method is a custom function
        aggre_func <- method
        method <- "custom"
    } else {
        method <- match.arg(method)
    }

    if (!is.matrix(geneList)) {
        stop("geneList should be a matrix of genes*models (true model or null models)")
    }

    aggre_func <- switch(method,
        mean = {
            function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                res <- mean(genelist[hits])
                return(res)
            }
        },
        median = {
            function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                res <- median(genelist[hits])
                return(res)
            }
        },
        meanabs = {
            function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                res <- mean(abs(genelist[hits]))
            }
        },
        meansqr = {
            function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                res <- mean(genelist[hits]^2)
                return(res)
            }
        },
        maxmean = {
            function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                X <- genelist[hits]
                # Mean of positive numbers
                pos.mean <- mean(X[X > 0])
                # Mean of negative numbers
                neg.mean <- mean(X[X < 0])
                if (is.na(pos.mean)) {
                    pos.mean <- 0
                } # when all value are neg
                if (is.na(neg.mean)) {
                    neg.mean <- 0
                } # when all value are pos
                res <- ifelse(pos.mean > abs(neg.mean), pos.mean, neg.mean)
                return(res)
            }
        },
        ks_orig = {
            function(genelist, geneSet) {
                # code obtained from DOSE: https://rdrr.io/bioc/DOSE/src/R/gsea.R
                genelist <- sort(genelist, decreasing = TRUE)
                geneSet <- intersect(geneSet, names(genelist))
                N <- length(genelist)
                Nh <- length(geneSet)
                Phit <- Pmiss <- numeric(N)
                hits <- names(genelist) %in% geneSet
                Phit[hits] <- abs(genelist[hits])^0 # raw rank
                NR <- sum(Phit)
                Phit <- cumsum(Phit / NR)
                Pmiss[!hits] <- 1 / (N - Nh)
                Pmiss <- cumsum(Pmiss)
                runningES <- Phit - Pmiss
                max.ES <- max(runningES)
                min.ES <- min(runningES)
                if (abs(max.ES) > abs(min.ES)) {
                    res <- max.ES
                } else {
                    res <- min.ES
                }
                return(res)
            }
        },
        ks_weighted = {
            function(genelist, geneSet) {
                # code obtained from DOSE: https://rdrr.io/bioc/DOSE/src/R/gsea.R
                genelist <- sort(genelist, decreasing = TRUE)
                geneSet <- intersect(geneSet, names(genelist))
                N <- length(genelist)
                Nh <- length(geneSet)
                Phit <- Pmiss <- numeric(N)
                hits <- names(genelist) %in% geneSet
                Phit[hits] <- abs(genelist[hits])^1 # weighted rank
                NR <- sum(Phit)
                Phit <- cumsum(Phit / NR)
                Pmiss[!hits] <- 1 / (N - Nh)
                Pmiss <- cumsum(Pmiss)
                runningES <- Phit - Pmiss
                max.ES <- max(runningES)
                min.ES <- min(runningES)
                if (abs(max.ES) > abs(min.ES)) {
                    res <- max.ES
                } else {
                    res <- min.ES
                }
                return(res)
            }
        },
        ks_pos_neg_sum = {
            function(genelist, geneSet) {
                # function borrowed from clusterprofiler
                genelist <- sort(genelist, decreasing = TRUE)
                geneSet <- intersect(geneSet, names(genelist))
                N <- length(genelist)
                Nh <- length(geneSet)
                Phit <- Pmiss <- numeric(N)
                hits <- names(genelist) %in% geneSet
                Phit[hits] <- abs(genelist[hits])^1 # weighted rank
                NR <- sum(Phit)
                Phit <- cumsum(Phit / NR)
                Pmiss[!hits] <- 1 / (N - Nh)
                Pmiss <- cumsum(Pmiss)
                runningES <- Phit - Pmiss
                max.ES <- max(runningES)
                min.ES <- min(runningES)
                res <- max.ES + min.ES
                return(res)
            }
        },
        sign_test = {
            function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                X <- genelist[hits]
                n.pos <- sum(X > 0)
                n.neg <- sum(X < 0)
                n.smaller <- ifelse(n.pos < n.neg, n.pos, n.neg)
                n.total <- n.pos + n.neg
                res <- ifelse(n.total <= 25, n.smaller, ((n.smaller + 0.5) - (n.total / 2)) / (sqrt(n.total) / 2))
                return(res)
            }
        },
        local_fdr = {
            function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet ## logical
                df <- data.frame(vals = as.numeric(genelist), hits = as.numeric(hits))
                fit.pos <- glm(hits ~ vals, family = binomial(link = "logit"), data = df[df$vals >= 0, ])
                fit.neg <- glm(hits ~ vals, family = binomial(link = "logit"), data = df[df$vals < 0, ])
                S.pos <- coef(summary(fit.pos))["vals", "z value"]
                S.neg <- coef(summary(fit.neg))["vals", "z value"]
                res <- ifelse(abs(S.pos) > abs(S.neg), S.pos, S.neg)
                return(res)
            }
        },
        rank_sum = {
            function(genelist, geneSet) {
                abslist.sorted <- sort(abs(genelist), decreasing = F) # sort by abs value
                ranklist <- c(1:length(genelist)) # create rank
                names(ranklist) <- names(abslist.sorted) # set the gene names
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                pos <- genelist > 0
                neg <- genelist < 0
                pos.name <- names(genelist)[hits & pos]
                neg.name <- names(genelist)[hits & neg]
                pos.ranksum <- sum(ranklist[pos.name])
                neg.ranksum <- sum(ranklist[neg.name])
                res <- ifelse(pos.ranksum < neg.ranksum, pos.ranksum, neg.ranksum)
                return(res)
            }
        },
        custom = {
            aggre_func
        },
        stop("Invalid aggregation method name")
    )


    if (!is.function(aggre_func)) {
        stop("Invalid method name")
    }

    gs_score <- sapply(1:ncol(geneList), function(x) {
        aggre_func(genelist = geneList[, x], geneSet = geneSet)
    })
    return(gs_score)
}



#' Aggregate Gene Set List
#'
#' This function aggregates multiple gene sets of interest after filtering, using a provided gene list.
#'
#' @param geneSetList A list of multiple gene sets of interest after filtering.
#' @param geneList A list of genes, which can be either a true gene list or a null gene list.
#' @param ... Additional arguments passed to the `aggregate_geneSet` function.
#'
#' @return A list of aggregated scores for each gene set.
#' @export
#'
#' @examples

aggregate_geneSetList <- function(geneList, geneSetList,  ...) {
  allgs.scores <- lapply(geneSetList, function(gs) {
    aggregate_geneSet(geneList = geneList, geneSet = gs, ...)
  })
  return(allgs.scores)
}



#' Aggregate Gene Set List Matching Co-Expression
#'
#' This function aggregates scores for multiple gene sets of interest, using a true gene list and sampled gene sets
#' while ensuring matching co-expression patterns.
#'
#' @param geneSetList A list of gene sets of interest after filtering.
#' @param sampled_geneSetList A list of sampled gene sets.
#' @param geneList.true A matrix representing the true gene list, with dimensions m x 1.
#' @param method method passed to `aggregate_geneSet`.
#' @return A list of aggregated scores for each gene set.
#' @examples
#' # Example usage with dummy data
#' geneList.true <- matrix(1:10, ncol = 1)
#' rownames(geneList.true) <- letters[1:10]
#' geneSetList <- list(set1 = c("a", "b", "c"), set2 = c("d", "e", "f"))
#' sampled_geneSetList <- list(set1 = list(c("g", "h", "i")), set2 = list(c("j", "k", "l")))
#' aggregate_geneSetList_matching_coexp(geneSetList, sampled_geneSetList, geneList.true)
#' @export
aggregate_geneSetList_matching_coexp <- function(geneList.true, 
                                                 geneSetList, 
                                                 sampled_geneSetList, 
                                                 method) {
  # Ensure geneList.true is a matrix with one column
  if (!is.matrix(geneList.true) || ncol(geneList.true) != 1) {
    stop('geneList.true should be a m x 1 matrix. Please include drop=FALSE when subsetting.')
  }
  
  # Ensure the order in geneSetList and sampled_geneSetList are the same
  if (!identical(names(geneSetList), names(sampled_geneSetList))) {
    stop('geneSetList and sampled_geneSetList are not matched.')
  }
  
  # Aggregate scores using mapply for parallel processing
  allgs.scores <- mapply(function(gs, sampled_gs) {
    geneList.null <- swap_geneList(geneList.true = geneList.true,
                                   orig_gs = gs,
                                   sampled_gs = sampled_gs)
    gs.score <- aggregate_geneSet(geneList = geneList.null,
                                  geneSet = gs,
                                  method=method)
    return(gs.score)
  }, geneSetList, sampled_geneSetList, SIMPLIFY = FALSE)
  
  return(allgs.scores)
}
