## Plotting layer (merged from the former IWTgg repository).

library(fastIWTlm)

load(system.file("extdata", "test_model.rdata", package = "fastIWTlm"))

## Draw all four panels (functional data + F-test, coefficients + t-tests).
autoplot(test_model)

## Capture the ggplot objects and compose them into side-by-side panels.
pl <- autoplot(test_model, save_output = TRUE, prompt = FALSE)
grids <- arrange_IWTlm(pl)
grids$data_global_plot
grids$coeff_global_plot
