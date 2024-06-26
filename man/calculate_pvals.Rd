% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{calculate_pvals}
\alias{calculate_pvals}
\title{Calculate P-Values}
\usage{
calculate_pvals(
  statList.true,
  statList.null,
  method = c("standard", "split_pos_neg")
)
}
\arguments{
\item{statList.true}{A named list of true statistics.}

\item{statList.null}{A named list of null statistics corresponding to the true statistics.}

\item{method}{The method to be used for p-value calculation. Either 'standard' or 'split_pos_neg'. Default is 'standard'.}
}
\value{
A list of calculated p-values.
}
\description{
This function calculates p-values based on the provided true statistics and null statistics lists.
It supports two methods for p-value calculation: 'standard' and 'split_pos_neg'.
}
