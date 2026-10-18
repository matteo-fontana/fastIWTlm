# fastIWTlm

**Fast Interval-Wise Testing (IWT) for functional regression — with a ggplot2 visualisation layer.**

This package merges two former repositories into one installable R package:

- **fastIWTlm** — the computation engine for functional regression IWT;
- **IWTgg** — the `ggplot2`-based plotting layer for `IWTlm` objects.

Now `autoplot()` works directly on the output of `fastIWTlm()`.

**One- and two-sample tests** (`iwt1()`, `iwt2()`, `twt2()`, `global2()`)
live in the companion package
[fdatestcpp](https://github.com/matteo-fontana/fdatestcpp) — see
[below](#one--and-two-sample-tests-use-fdatestcpp).

## Features

| Function | Purpose | Backend |
|---|---|---|
| `fastIWTlm(formula, ...)` | IWT for functional-on-scalar / functional-on-functional regression | R, or **RcppEigen + OpenMP** via `engine = "rcpp"` |
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
```

## One- and two-sample tests: use fdatestcpp

Earlier versions of this package shipped a two-sample engine (`fastIWT2()`,
with the pure-R reference `IWT2()` and the `IWTimage()` heatmap). They have
been removed. Use [fdatestcpp](https://github.com/matteo-fontana/fdatestcpp)
instead: a C++ port of the development version of
[fdatest](https://github.com/permaverse/fdatest) with the same interface.
It is reproducible with `set.seed()`, gives results `identical()` to fdatest,
and is much faster.

```r
# remotes::install_github("matteo-fontana/fdatestcpp")
library(fdatestcpp)

# data1, data2: one curve per row, one grid point per column
res <- iwt2(data1, data2, mu = 0, n_perm = 2000, recycle = FALSE)
res$adjusted_pvalues
plot(res)
```

Mapping from the removed functions:

| Removed | Use instead |
|---|---|
| `fastIWT2(Data1, Data2, Mu, B, alt, maxrow = 0, paired, recycle)` | `fdatestcpp::iwt2(data1, data2, mu, n_perm = B, alternative = alt, paired = paired, recycle = recycle)` |
| `IWT2(data1, data2, mu, B, maxrow = 1, paired, recycle, alternative)` | the same `fdatestcpp::iwt2()` call |
| `IWTimage(result)` | `plot(res)` / `autoplot(res)` on the `iwt2()` result |

`fdatestcpp::iwt2()` has no `maxrow` argument: it always computes the full
set of intervals, as `maxrow = 0` (`fastIWT2`) or `maxrow = 1` (`IWT2`) did.

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
R/     package R sources (regression, plotting, Rcpp glue)
src/   C++/Rcpp engine: fastIWTlmperm.cpp (+ Makevars)
man/   documentation
inst/extdata/       example datasets for fastIWTlm()
inst/examples/      runnable demo scripts
docs/              abstract, report and presentation (not built)
```

## Status

Built, installed and checked with **R 4.3.3** (RcppEigen 0.3.4, Rcpp 1.0.12,
ggplot2 3.4.4, data.table 1.14.10):

- `R CMD check`: **0 errors, 0 warnings**, only routine NOTEs (compiled-object
  size, `GNU make` / OpenMP `SystemRequirements`, and `ggpubr` when that
  optional Suggests is absent).
- Re-checked after removing the two-sample engine, with
  `R CMD check --as-cran` on R 4.6.1 (Windows): 0 errors, 0 warnings,
  3 NOTEs (new submission, `pandoc` not installed locally, and a stray `NULL`
  directory in the check folder that no example creates, since all are
  `\dontrun`). The `rand()` and `class()` comparison NOTEs are gone.
- `Rcpp::compileAttributes()` reproduces `src/RcppExports.cpp` verbatim.
- **Engine parity:** under a common `set.seed()` the `"R"` and `"rcpp"`
  engines are numerically identical (max abs difference `0` on the bundled
  example), with the RcppEigen kernel ~25–30x faster on that example.

Reproduce:

```r
Rcpp::compileAttributes()          # no diff in RcppExports.*
# R CMD build . && R CMD check fastIWTlm_0.1.0.tar.gz
source(system.file("examples", "fastIWTlm_test_code.R", package = "fastIWTlm"))
```
