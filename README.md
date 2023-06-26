
<!-- README.md is generated from README.Rmd. Please edit that file -->

# segregatr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/segregatr)](https://CRAN.R-project.org/package=segregatr)
[![](https://cranlogs.r-pkg.org/badges/grand-total/segregatr?color=yellow)](https://cran.r-project.org/package=segregatr)
[![](https://cranlogs.r-pkg.org/badges/last-month/segregatr?color=yellow)](https://cran.r-project.org/package=segregatr)
<!-- badges: end -->

The goal of **segregatr** is to provide segregation analysis for
clinical variant classification.

## Installation

You can install **segregatr** from CRAN as follows:

``` r
install.packages("segregatr")
```

Alternatively, obtain the latest development version from GitHub:

``` r
devtools::install_github("magnusdv/segregatr")
```

## Example

``` r
library(segregatr)
```

The family below shows four brothers, all affected with a rare dominant
disease with 90% penetrance and phenocopy rate 1%. The parents have
unknown affection status. All four brothers are shown to carry a
candidate variant.

<img src="man/figures/README-sibex-1.png" width="40%" style="display: block; margin: auto;" />

We will use **segregatr** to analyse the co-segregation of the variant
and the disease in this pedigree. Specifically we want to compute the
*full-likelihood Bayes factor* (FLB), quantifying the evidence that the
variant is pathogenic.

To create the pedigree we use the `nuclearPed()` function from the
**pedtools** package, which is automatically loaded together with
**segregatr**.

``` r
x = nuclearPed(4)
```

Then we run the `FLB()` function, filling in the necessary data:

``` r
FLB(x, carriers = 3:6, affected = 3:6, unknown = 1:2,
    freq = 0.0001, penetrances = c(0.01, 0.9, 0.9), proband = 3)
#> [1] 7.732161
```

The resulting FLB score is less than 8, which unfortunately only
indicates suggestive evidence for pathogenicity.
