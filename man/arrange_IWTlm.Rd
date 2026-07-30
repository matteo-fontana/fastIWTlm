\name{arrange_IWTlm}
\alias{arrange_IWTlm}
\title{Compose IWTlm ggplot Panels}
\description{
  Arranges the four \code{ggplot} objects returned by
  \code{\link{autoplot.IWTlm}} into two side-by-side composite panels:
  functional data next to the F-test p-value function, and functional
  coefficients next to the t-test p-value functions.
}
\usage{
arrange_IWTlm(plot_list)
}
\arguments{
  \item{plot_list}{the list returned by
    \code{autoplot.IWTlm(..., save_output = TRUE)}, containing
    \code{data_plot}, \code{coeff_plot}, \code{pval_t_plot} and
    \code{pval_plot_F}.}
}
\value{
  A list with two \code{ggarrange} objects: \code{data_global_plot} and
  \code{coeff_global_plot}.
}
\details{Requires the \pkg{ggpubr} package.}
\seealso{\code{\link{autoplot.IWTlm}}}
