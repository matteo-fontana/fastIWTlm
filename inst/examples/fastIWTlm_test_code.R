## Functional-on-scalar/functional regression via Interval-Wise Testing,
## with a parity check between the pure-R and the Rcpp permutation engines.

library(fastIWTlm)

load(system.file("extdata", "test_data_IWTlm.rdata", package = "fastIWTlm"))
attach(cov)

f <- co2_mat_log ~ population_mat_log + gdppc_mat_log + end + ff + lc

## Headline (pure-R) engine
test_model <- fastIWTlm(f, B = 1000)

## Same procedure on the parallel RcppEigen engine. Under a common seed the two
## engines draw the same permutation matrix and must agree to FP tolerance.
set.seed(1); m_R   <- fastIWTlm(f, B = 500, engine = "R")
set.seed(1); m_cpp <- fastIWTlm(f, B = 500, engine = "rcpp", threads = 4)
stopifnot(isTRUE(all.equal(m_R$adjusted_pval_part, m_cpp$adjusted_pval_part,
                           tol = 1e-8)))

## Visualise
autoplot(test_model)
