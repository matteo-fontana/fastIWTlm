# fastIWTlm

**Fast Interval-Wise Testing (IWT) for functional linear models.**

fastIWTlm tests linear models with a functional response and scalar or
functional (concurrent) covariates. It returns the global F test and the
partial t tests as functions over the domain, with p-values adjusted to
control the interval-wise error rate (Pini and Vantini, 2017; Abramowicz et
al., 2018), and plots the results with ggplot2.

## Origin of the code

`fastIWTlm()` is adapted from R code by Alessia Pini for interval-wise
testing in linear models with functional covariates. That code builds on the
interval-wise testing code for functional-on-scalar linear models by Alessia
Pini and Simone Vantini (`ITPlmbspline()` in the
[fdatest](https://CRAN.R-project.org/package=fdatest) package). Matteo
Fontana packaged it and added the RcppEigen permutation kernel and the
ggplot2 plots.

## Features

| Function | Purpose | Backend |
|---|---|---|
| `fastIWTlm(formula, ...)` | IWT for functional-on-scalar / functional-on-functional regression | R, or **RcppEigen + OpenMP** via `engine = "rcpp"` |
| `autoplot(model)` / `plot(model)` | Four-panel diagnostic plots | ggplot2 |
| `arrange_IWTlm(plot_list)` | Compose the panels side by side | ggplot2 / ggpubr |
| `summary(model)` | Tabular summary | — |

## Installation

The package is compiled from source: on Windows this needs
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), on macOS the Xcode
command line tools. OpenMP is used when the compiler supports it.

```r
# install.packages("remotes")
remotes::install_github("matteo-fontana/fastIWTlm")
```

## Quick start

```r
library(fastIWTlm)

load(system.file("extdata", "test_data_IWTlm.rdata", package = "fastIWTlm"))
attach(cov)
f <- co2_mat_log ~ population_mat_log + gdppc_mat_log + end + ff + lc

model <- fastIWTlm(f, B = 1000)          # pure-R engine
summary(model)
autoplot(model)                          # ggplot2 diagnostics

## Accelerate the permutation loop (parallel RcppEigen kernel):
model2 <- fastIWTlm(f, B = 1000, engine = "rcpp", threads = 2)
```

## The `engine = "rcpp"` kernel

The permutation loop refits the model for every partial test at every
abscissa point on each of `B` permutations. Because the per-point design
matrix is **fixed across permutations** (only the response is permuted or
reconstructed), the projection `(X'X)^{-1}X'` and `diag((X'X)^{-1})` can be
computed once and reused. `engine = "rcpp"` does exactly this in a parallel
RcppEigen kernel (`src/fastIWTlmperm.cpp`).

The permutation indices are drawn in R and passed to the kernel, so under a
common `set.seed()` the two engines give the same results:

```r
set.seed(1); mR   <- fastIWTlm(f, B = 500, engine = "R")
set.seed(1); mCpp <- fastIWTlm(f, B = 500, engine = "rcpp", threads = 2)
all.equal(mR$adjusted_pval_part, mCpp$adjusted_pval_part)   # TRUE
```

## Repository layout

```
R/              package R sources (regression, plotting, Rcpp glue)
src/            C++/Rcpp kernel: fastIWTlmperm.cpp (+ Makevars)
man/            documentation
inst/extdata/   example datasets
inst/examples/  runnable demo scripts
```

## Authors

- Matteo Fontana (maintainer)
- Alessia Pini — R code from which `fastIWTlm()` is adapted
- Simone Vantini — co-author of the fdatest code that code builds on
- Aymeric Stamm — contributor; maintainer of
  [fdatest](https://github.com/permaverse/fdatest)

## References

Pini, A. and Vantini, S. (2017). Interval-wise testing for functional data.
*Journal of Nonparametric Statistics*, 29(2), 407–424.
<https://doi.org/10.1080/10485252.2017.1306627>

Abramowicz, K., Häger, C. K., Pini, A., Schelin, L., Sjöstedt de Luna, S. and
Vantini, S. (2018). Nonparametric inference for functional-on-scalar linear
models applied to knee kinematic hop data after injury of the anterior
cruciate ligament. *Scandinavian Journal of Statistics*, 45(4), 1036–1061.
<https://doi.org/10.1111/sjos.12333>

## License

GPL-3
