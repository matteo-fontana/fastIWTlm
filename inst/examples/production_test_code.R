## Two-sample Interval-Wise Testing with the C++/Rcpp engine (fastIWT2).

library(fastIWTlm)

load(system.file("extdata", "test-data.rdata", package = "fastIWTlm"))

## Directly call the compiled, OpenMP-parallelised engine.
results <- fastIWT2(data1, data2, mu, B = 5000, MPthreads = 4)

str(results)
