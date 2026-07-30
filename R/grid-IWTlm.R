# Compose the four ggplot outputs of autoplot.IWTlm() into two side-by-side
# panels: (functional data | F-test p-values) and (coefficients | t-test
# p-values). `plot_list` is the list returned by
# autoplot.IWTlm(..., save_output = TRUE).
arrange_IWTlm <- function(plot_list) {
  if (!requireNamespace("ggpubr", quietly = TRUE))
    stop("Package 'ggpubr' is required for arrange_IWTlm(); please install it.")

  data_global_plot <- ggpubr::ggarrange(
    plot_list$data_plot, plot_list$pval_plot_F,
    ncol = 2, legend = 'bottom')

  coeff_global_plot <- ggpubr::ggarrange(
    plot_list$coeff_plot, plot_list$pval_t_plot,
    ncol = 2, common.legend = TRUE, legend = 'bottom')

  list(data_global_plot  = data_global_plot,
       coeff_global_plot = coeff_global_plot)
}
