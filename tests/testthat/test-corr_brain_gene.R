library(testthat)
library(BrainEnrich)

# Set seed for reproducibility
set.seed(123)

# Number of samples and genes/regions
n_samples <- 1
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


test_that("corr_brain_gene works with Pearson correlation", {
  result <- corr_brain_gene(gene_data, brain_data, method = "pearson")
  # Check if the result is a matrix and its dimensions match the brain_data
  expect_true(is.matrix(result), "Result should be a matrix")
  expect_equal(nrow(result), ncol(gene_data)) # result should be gene x sample
  expect_equal(ncol(result), ncol(brain_data)) # result should be gene x sample
  expect_equal(rownames(result), colnames(gene_data))
  expect_equal(colnames(result), colnames(brain_data))
})

test_that("corr_brain_gene works with Spearman correlation", {
  result <- corr_brain_gene(gene_data, brain_data, method = "spearman")
  # Check if the result is a matrix and its dimensions match the brain_data
  expect_true(is.matrix(result), "Result should be a matrix")
  expect_equal(nrow(result), ncol(gene_data)) # result should be gene x sample
  expect_equal(ncol(result), ncol(brain_data)) # result should be gene x sample
  expect_equal(rownames(result), colnames(gene_data))
  expect_equal(colnames(result), colnames(brain_data))
})

test_that("corr_brain_gene works with PLS1 component regression", {
  result <- corr_brain_gene(gene_data, brain_data, method = "pls1c")
  expect_true(is.matrix(result), "Result should be a matrix")
  expect_equal(nrow(result), ncol(gene_data)) # result should be gene x sample
  expect_equal(ncol(result), ncol(brain_data)) # result should be gene x sample
  expect_equal(rownames(result), colnames(gene_data))
  expect_equal(colnames(result), colnames(brain_data))
})

test_that("corr_brain_gene works with PLS1 weights", {
  result <- corr_brain_gene(gene_data, brain_data, method = "pls1w")
  expect_true(is.matrix(result), "Result should be a matrix")
  expect_equal(nrow(result), ncol(gene_data)) # result should be gene x sample
  expect_equal(ncol(result), ncol(brain_data)) # result should be gene x sample
  expect_equal(rownames(result), colnames(gene_data))
  expect_equal(colnames(result), colnames(brain_data))
})

test_that("corr_brain_gene works with custom correlation function", {
  custom_function <- function(gene, brain) {
    abs(cor(gene, brain, method = "pearson") - 0.1)
  }
  result <- corr_brain_gene(gene_data, brain_data, method = custom_function)
  expect_true(is.matrix(result), "Result should be a matrix")
  expect_equal(nrow(result), ncol(gene_data)) # result should be gene x sample
  expect_equal(ncol(result), ncol(brain_data)) # result should be gene x sample
  expect_equal(rownames(result), colnames(gene_data))
  expect_equal(colnames(result), colnames(brain_data))

  result1 <- corr_brain_gene(gene_data, brain_data, method = "pearson", r2z = FALSE)
  result2 <- abs(result1 - 0.1)
  # expect_equivalent(result1, result3, tolerance = 1e-10)
  expect_equal(result, result2, tolerance = 1e-10, ignore_attr = TRUE) # attribute is_fisherz is not the same
})
