\name{summary.IWTlm}
\alias{summary.IWTlm}
\title{Summarise an IWTlm Object}
\description{
  Tabular summary of an \code{"IWTlm"} object: the minimum interval-wise
  adjusted p-value of each partial (t) test with significance codes, the
  range of the functional R-squared, and the minimum adjusted p-value of
  the global (F) test.
}
\usage{
\method{summary}{IWTlm}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{"IWTlm"} returned by
    \code{\link{fastIWTlm}}.}
  \item{\dots}{currently ignored.}
}
\value{
  A list with elements \code{call}, \code{ttest}, \code{R2} and
  \code{ftest}.
}
\seealso{\code{\link{fastIWTlm}}}
