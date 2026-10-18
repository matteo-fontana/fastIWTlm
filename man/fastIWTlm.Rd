\name{fastIWTlm}
\alias{fastIWTlm}
\title{Interval-Wise Testing for Functional Regression}
\description{
  Interval-Wise Testing (IWT) for functional-on-scalar and
  functional-on-functional linear regression models. A pointwise
  permutation scheme produces global (F-type) and partial (t-type)
  test-statistic functions, which are then combined over intervals and
  adjusted, controlling the interval-wise error rate.

  The permutation loop can optionally be dispatched to a parallel
  RcppEigen kernel via \code{engine = "rcpp"} for a large speed-up on the
  dominant cost of the procedure.
}
\usage{
fastIWTlm(formula, B = 1000, method = "residuals",
          engine = c("R", "rcpp"), threads = 1)
}
\arguments{
  \item{formula}{a model \code{\link{formula}}. The response must be a
    functional dataset (an \eqn{n \times p} matrix); covariates may be
    functional (\eqn{n \times p} matrices) or scalar (length-\eqn{n}
    vectors), e.g. \code{y ~ x1 + x2}.}
  \item{B}{number of permutations.}
  \item{method}{permutation scheme for the partial tests, either
    \code{"residuals"} (default) or \code{"responses"}.}
  \item{engine}{\code{"R"} (default) runs the pure-R permutation loop;
    \code{"rcpp"} dispatches it to a parallel RcppEigen kernel. The two
    engines are numerically equivalent under a common \code{\link{set.seed}}.
    \code{"rcpp"} requires at least one covariate (\code{nvar > 0}); for
    intercept-only models it falls back to \code{"R"} with a warning.}
  \item{threads}{number of OpenMP threads used when \code{engine = "rcpp"}.}
}
\value{
  An object of class \code{"IWTlm"}: a list with, among others,
  \code{unadjusted_pval_F}, \code{adjusted_pval_F},
  \code{unadjusted_pval_part}, \code{adjusted_pval_part},
  \code{coeff.regr.eval}, \code{fitted.eval}, \code{residuals.eval},
  \code{R2.eval} and the heatmap matrices. Use \code{\link{summary.IWTlm}}
  and \code{\link{autoplot.IWTlm}} to inspect it.
}
\seealso{\code{\link{summary.IWTlm}}, \code{\link{autoplot.IWTlm}},
  \code{\link{arrange_IWTlm}}; \pkg{fdatestcpp}
  (\url{https://github.com/matteo-fontana/fdatestcpp}) for one- and
  two-sample interval-wise tests.}
\examples{
\dontrun{
load(system.file("extdata", "test_data_IWTlm.rdata", package = "fastIWTlm"))
attach(cov)
f <- co2_mat_log ~ population_mat_log + gdppc_mat_log + end + ff + lc

set.seed(1); m_R   <- fastIWTlm(f, B = 500, engine = "R")
set.seed(1); m_cpp <- fastIWTlm(f, B = 500, engine = "rcpp", threads = 4)
all.equal(m_R$adjusted_pval_part, m_cpp$adjusted_pval_part)

autoplot(m_R)
}
}
