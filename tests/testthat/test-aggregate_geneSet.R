library(testthat)
library(BrainEnrich)

# Set seed for reproducibility
set.seed(123)

# Number of samples and genes/regions
n_samples <- 3
n_genes <- 55
n_regions <- 33

# Simulate gene expression data
gene_data <- mvtnorm::rmvnorm(n_regions, mean = rep(0, n_genes))
rownames(gene_data) <- paste("Region", 1:n_regions, sep = "_")
colnames(gene_data) <- paste("Gene", 1:n_genes, sep = "_")

# Simulate brain activity data with some correlation to gene data
brain_data <- mvtnorm::rmvnorm(n_regions, mean = rep(0, n_samples))
rownames(brain_data) <- paste("Region", 1:n_regions, sep = "_")
colnames(brain_data) <- paste("Sample", 1:n_samples, sep = "_")
geneList=corr_brain_gene(gene_data, brain_data, method = "pearson")

# simulate a vector of gene names
geneSet <- paste("Gene", sample(1:10, 5), sep="_")


test_that("aggregate_geneSet works with mean method", {
  result <- aggregate_geneSet(geneList, geneSet, method = "mean")
  expected_result <- apply(geneList[geneSet, , drop = FALSE],2,mean)
  expect_equal(result, as.numeric(expected_result), tolerance = 1e-8)
})

test_that("aggregate_geneSet works with median method", {
  result <- aggregate_geneSet(geneList, geneSet, method = "median")
  expected_result <- apply(geneList[geneSet, , drop = FALSE],2,median)
  expect_equal(result, as.numeric(expected_result), tolerance = 1e-8)
})


test_that("aggregate_geneSet works with meanabs method", {
  result <- aggregate_geneSet(geneList, geneSet, method = "meanabs")
  expected_result <- apply(geneList[geneSet, , drop = FALSE],2,function(x) mean(abs(x)))
  expect_equal(result, as.numeric(expected_result), tolerance = 1e-8)
})

test_that("aggregate_geneSet handles custom function", {
  custom_func <- function(genelist, geneSet) {
                geneSet <- intersect(geneSet, names(genelist))
                hits <- names(genelist) %in% geneSet
                res <- mean(genelist[hits])+ 10
                return(res)
  }
  result <- aggregate_geneSet(geneList, geneSet, method = custom_func)
  expected_result <- aggregate_geneSet(geneList, geneSet, method = 'mean') + 10
  expect_equal(result, expected_result, tolerance = 1e-8)
})

