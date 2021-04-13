
<!-- README.md is generated from README.Rmd. Please edit that file -->

# segregatr

<!-- badges: start -->
<!-- badges: end -->

The goal of **segregatr** is to provide segregation analysis for
clinical variant classification.

## Installation

You can install the development version of **segregatr** from GitHub
with:

``` r
# install.packages("devtools")
devtools::install_github("magnusdv/segregatr")
```

``` r
library(segregatr)
#> Loading required package: pedtools
```

## Example

The family below shows four brothers, all affected with a rare dominant
disease with 90% penetrance and phenocopy rate 1%. The parents have
unknown affection status. All four brothers are shown to carry a
candidate variant, warranting a segregation analysis. pathogenic
variant.

<img src="man/figures/README-sibex-1.png" width="40%" />

In order to compute the full-likelihood Bayes factor, we first create
the pedigree.

``` r
x = nuclearPed(4)
```

Then we run the `FLB()` function, filling in the necessary data:

``` r
FLB(x, carriers = 3:6, aff = 3:6, unknown = 1:2,
    freq = 0.0001, penetrances = c(0.01, 0.9, 0.9), proband = 3)
#> [1] 7.732161
```

The answer indicates only suggestive evidence for pathogenicity.
