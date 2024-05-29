# install necessary packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}

if (!requireNamespace("roxygen2", quietly = TRUE)) {
  install.packages("roxygen2")
}

if (!requireNamespace("testthat", quietly = TRUE)) {
  install.packages("testthat")
}

devtools::has_devel()

devtools::document()
devtools::build()
devtools::install()
usethis::use_mit_license()
usethis::use_citation()
usethis::use_agpl3_license()



usethis::use_data()

usethis::use_testthat() # Create testthat infrastructure
usethis::use_test('corr_brain_gene') # Create test file for corr_brain_gene function
devtools::test(filter = "corr_brain_gene") # Run tests for corr_brain_gene function
# testthat::test_file("tests/testthat/test-corr_brain_gene.R") # Run tests for corr_brain_gene function

usethis::use_test('aggregate_geneSet') # Create test file for aggregate_gene_data function
devtools::test(filter = "aggregate_geneSet") 

usethis::use_test('get_geneSets.R')
usethis::use_test('brainenrich.R')
usethis::use_test('brainscore.R')

usethis::use_vignette('brainenrich')
usethis::use_vignette('visualization')
usethis::use_vignette('brainscore')

usethis::use_pkgdown()
usethis::edit_r_environ()
usethis::browse_github_pat()

usethis::use_pkgdown_github_pages()
pkgdown::build_site()
# prepare example data

