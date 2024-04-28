library(MASS)  # For simulating bimodal distribution
library(BrainTransR)
set.seed(123)  # Set seed for reproducibility



# Number of samples and genes/regions
n_samples <- 1
n_genes <- 550
n_regions <- 33

# Simulate gene expression data
gene_data <- mvtnorm::rmvnorm(n_regions, mean = rep(0, n_genes))
rownames(gene_data) <- paste("Region", 1:n_regions, sep = "_")
colnames(gene_data) <- paste("Gene", 1:n_genes, sep = "_")

# Simulate brain activity data with some correlation to gene data
brain_data <- mvtnorm::rmvnorm(n_regions, mean = rep(0, n_samples))
rownames(brain_data) <- paste("Region", 1:n_regions, sep = "_")
colnames(brain_data) <- paste("Sample", 1:n_samples, sep = "_")


result <- cor(gene_data, brain_data, method = "pearson")
d=diptest::dip.test(result)
