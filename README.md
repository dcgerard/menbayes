
<!-- README.md is generated from README.Rmd. Please edit that file -->

# menbayes

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/660268756.svg)](https://zenodo.org/doi/10.5281/zenodo.12189055)
[![NSF-2132247](https://img.shields.io/badge/NSF-2132247-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2132247)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/dcgerard/menbayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dcgerard/menbayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Provides a suite of tests for segregation distortion in F1 polyploid
populations (for now, just tetraploids). This is under different
assumptions of meiosis. The main functions are:

- `bayes_men_g4()`: Bayesian test for segregation distortion using known
  genotypes.
- `bayes_men_gl4()`: Bayesian test for segregation distortion using
  genotype likelihoods.
- `lrt_men_g4()`: Likelihood ratio test for segregation distortion using
  known genotypes.
- `lrt_men_gl4()`: Likelihood ratio test for segregation distortion
  using genotype likelihoods.
- `offspring_gf_2()`: Offspring genotype frequencies under the two
  parameter model of meiosis.
- `offspring_gf_3()`: Offspring genotype frequencies under the three
  parameter model of meiosis.
- `simf1g()`: Simulate genotypes from an F1 population of tetraploids.
- `simf1gl()`: Simulate genotype likelihoods from an F1 population of
  tetraploids.

Details of these methods may be found in Gerard et al. (2025).

## Installation

You can install the development version of menbayes from
[GitHub](https://github.com/dcgerard/menbayes) with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/menbayes")
```

## Code of Conduct

Please note that the menbayes project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## References

Gerard D, Thakkar M, & Ferrão LFV (2025). “Tests for segregation
distortion in tetraploid F1 populations.” *Theoretical and Applied
Genetics*, *138*(30), p. 1–13.
[doi:10.1007/s00122-025-04816-z](https://doi.org/10.1007/s00122-025-04816-z).

## Acknowledgments

This material is based upon work supported by the National Science
Foundation under Grant No. 2132247.
