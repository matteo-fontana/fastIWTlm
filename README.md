# fastIWTlm

**Fast Interval-Wise Testing (IWT) for functional data — with a ggplot2 visualisation layer.**

This package merges two former repositories into one installable R package:

- **fastIWTlm** — the computation engine (functional regression IWT and a
  high-performance C++/Rcpp two-sample engine);
- **IWTgg** — the `ggplot2`-based plotting layer for `IWTlm` objects.

Now `autoplot()` works directly on the output of `fastIWTlm()`.

## Features

| Function | Purpose | Backend |
|---|---|---|
| `fastIWTlm(formula, ...)` | IWT for functional-on-scalar / functional-on-functional regression | R, or **RcppEigen + OpenMP** via `engine = "rcpp"` |
| `fastIWT2(Data1, Data2, Mu, ...)` | Two-sample IWT | **C++/Rcpp (RcppEigen + OpenMP)** |
| `IWT2(data1, data2, ...)` | Two-sample IWT | pure R (reference) |
| `IWTimage(result)` | p-value heatmap | base graphics |
| `autoplot(model)` / `plot(model)` | Four-panel diagnostic plots | ggplot2 |
| `arrange_IWTlm(plot_list)` | Compose the panels side by side | ggplot2 / ggpubr |
| `summary(model)` | Tabular summary | — |

## Installation

Requires a C++11 toolchain with OpenMP and the `Rcpp`, `RcppEigen`,
`ggplot2` and `data.table` packages.

```r
# install.packages(c("Rcpp", "RcppEigen", "ggplot2", "data.table"))
# from a local clone:
# R CMD INSTALL fastIWTlm
# or, with devtools:
devtools::install()
```

## Quick start

```r
library(fastIWTlm)

## --- Functional regression -------------------------------------------------
load(system.file("extdata", "test_data_IWTlm.rdata", package = "fastIWTlm"))
attach(cov)
f <- co2_mat_log ~ population_mat_log + gdppc_mat_log + end + ff + lc

model <- fastIWTlm(f, B = 1000)          # pure-R engine
summary(model)
autoplot(model)                          # ggplot2 diagnostics

## Accelerate the permutation loop (parallel RcppEigen kernel):
model2 <- fastIWTlm(f, B = 1000, engine = "rcpp", threads = 4)

## --- Two-sample test (compiled engine) ------------------------------------
load(system.file("extdata", "test-data.rdata", package = "fastIWTlm"))
res <- fastIWT2(data1, data2, mu, B = 5000, MPthreads = 4)
```

## The `engine = "rcpp"` extension

The regression function `fastIWTlm()` was historically pure R: its
permutation loop re-fits `lm.fit` for every partial test at every abscissa
point on each of `B` iterations. Because the per-point design matrix is
**fixed across permutations** (only the response is permuted or
reconstructed), the projection `(X'X)^{-1}X'` and `diag((X'X)^{-1})` can be
computed once and reused. `engine = "rcpp"` does exactly this in a
parallel RcppEigen kernel (`src/fastIWTlmperm.cpp`).

The permutation indices are drawn in R and passed to the kernel, so under a
common `set.seed()` the two engines are numerically identical:

```r
set.seed(1); mR   <- fastIWTlm(f, B = 500, engine = "R")
set.seed(1); mCpp <- fastIWTlm(f, B = 500, engine = "rcpp", threads = 4)
all.equal(mR$adjusted_pval_part, mCpp$adjusted_pval_part)   # TRUE
```

## Repository layout

```
R/     package R sources (regression, two-sample, plotting, heatmap, Rcpp glue)
src/   C++/Rcpp engines: fastIWT2.cpp, fastIWTlmperm.cpp (+ Makevars)
man/   documentation
inst/extdata/       example datasets and text fixtures
inst/examples/      runnable demo / benchmark scripts
inst/standalone-cpp/ the original non-Rcpp GetPot/Eigen command-line tool
dev/cpp/            experimental C++ kept for reference (not built)
docs/              abstract, report and presentation (not built)
```

## Notes for maintainers

The package was assembled in an environment without R, so the Rcpp glue
(`R/RcppExports.R`, `src/RcppExports.cpp`), `NAMESPACE` and `man/*.Rd` were
hand-authored to match the standard `Rcpp::compileAttributes()` / roxygen
output. Before release, run in an R session:

```r
Rcpp::compileAttributes()   # should reproduce RcppExports.* with no diff
# then:
# R CMD build . && R CMD check fastIWTlm_0.1.0.tar.gz
```
