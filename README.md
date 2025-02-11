
<!-- README.md is generated from README.Rmd. Please edit that file -->

# statgenIBD <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![](https://www.r-pkg.org/badges/version/statgenIBD)](https://www.r-pkg.org/pkg/statgenIBD)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/statgenIBD)](https://www.r-pkg.org/pkg/statgenIBD)
[![R-CMD-check](https://github.com/Biometris/statgenIBD/workflows/R-CMD-check/badge.svg)](https://github.com/Biometris/statgenIBD/actions?workflow=R-CMD-check)
[![codecov](https://codecov.io/gh/Biometris/statgenIBD/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Biometris/statgenIBD)

**statgenIBD** is an R package for calculating Identity By Descent (IBD)
probabilities for biparental, three and four-way crosses. Calculations
are based on Hidden Markov Models (HMM) and inheritance vectors. The HMM
calculations are implemented in `Rcpp/C++`.

For more complicated pedigrees
[RABBIT](https://github.com/Biometris/RABBIT) can be used to calculate
the IBDs.

## Installation

- Install from CRAN:

``` r
install.packages("statgenIBD")
```

- Install latest development version from GitHub (requires
  [remotes](https://github.com/r-lib/remotes) package):

``` r
remotes::install_github("Biometris/statgenIBD", ref = "develop", dependencies = TRUE)
```
