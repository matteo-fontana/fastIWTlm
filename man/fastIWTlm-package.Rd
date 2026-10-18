\name{fastIWTlm-package}
\alias{fastIWTlm-package}
\docType{package}
\title{Fast Interval-Wise Testing for Functional Linear Models}
\description{
  Interval-Wise Testing (IWT) for linear models with a functional response
  and scalar or functional (concurrent) covariates.
}
\details{
  Main entry points:
  \describe{
    \item{\code{\link{fastIWTlm}}}{IWT for functional linear models;
      returns an object of class \code{"IWTlm"}. Its permutation loop can
      be accelerated with \code{engine = "rcpp"} (RcppEigen + OpenMP).}
    \item{\code{\link{summary.IWTlm}}}{tabular summary of a fit.}
    \item{\code{\link{autoplot.IWTlm}}, \code{\link{arrange_IWTlm}}}{the
      \pkg{ggplot2} visualisation layer for \code{"IWTlm"} objects.}
  }

  Origin of the code: \code{\link{fastIWTlm}} is adapted from R code by
  Alessia Pini for interval-wise testing in linear models with functional
  covariates. That code builds on the interval-wise testing code for
  functional-on-scalar linear models by Alessia Pini and Simone Vantini
  (\code{ITPlmbspline()} in the \pkg{fdatest} package). Matteo Fontana
  packaged it and added the RcppEigen permutation kernel and the
  \pkg{ggplot2} plots.
}
\references{
  Pini, A. and Vantini, S. (2017). Interval-wise testing for functional
  data. \emph{Journal of Nonparametric Statistics}, 29(2), 407--424.
  \doi{10.1080/10485252.2017.1306627}

  Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S.
  and Vantini, S. (2018). Nonparametric inference for functional-on-scalar
  linear models applied to knee kinematic hop data after injury of the
  anterior cruciate ligament. \emph{Scandinavian Journal of Statistics},
  45(4), 1036--1061. \doi{10.1111/sjos.12333}
}
\author{
  Matteo Fontana (maintainer), Alessia Pini, Simone Vantini, and Aymeric
  Stamm (contributor; maintainer of \pkg{fdatest}). See Details for the
  origin of the code.
}
\keyword{package}
