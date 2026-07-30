# Declare non-standard-evaluation column names used by the data.table / ggplot2
# plotting layer so that R CMD check does not emit
# "no visible binding for global variable" NOTEs.
utils::globalVariables(c(
  "min_rect", "max_rect", "value", "significance",
  "variable", "abscissa_smooth", "type"
))
