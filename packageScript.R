library(usethis)
library(devtools)
library(roxygen2)

has_devel()

#usethis::create_package("BrainTransR")


devtools::document()
devtools::build()
usethis::use_mit_license()
