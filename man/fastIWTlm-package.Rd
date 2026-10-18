\name{fastIWTlm-package}
\alias{fastIWTlm-package}
\docType{package}
\title{Fast Interval-Wise Testing with ggplot2 Visualisation}
\description{
  \pkg{fastIWTlm} merges the former \code{fastIWTlm} (computation) and
  \code{IWTgg} (plotting) projects into a single package for Interval-Wise
  Testing (IWT) of functional data.
}
\details{
  Main entry points:
  \describe{
    \item{\code{\link{fastIWTlm}}}{IWT for functional regression models;
      returns an object of class \code{"IWTlm"}. Its permutation loop can
      be accelerated with \code{engine = "rcpp"} (RcppEigen + OpenMP).}
    \item{\code{\link{autoplot.IWTlm}}, \code{\link{arrange_IWTlm}}}{the
      \pkg{ggplot2} visualisation layer for \code{"IWTlm"} objects.}
  }
  One- and two-sample interval-wise tests are provided by the companion
  package \pkg{fdatestcpp}
  (\url{https://github.com/matteo-fontana/fdatestcpp}).
}
\author{Matteo Fontana and contributors.}
\keyword{package}
