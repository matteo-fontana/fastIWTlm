## Two-sample IWT: compare the pure-R and the compiled engines on the same data.

library(fastIWTlm)

ex <- function(f) system.file("extdata", f, package = "fastIWTlm")
data1 <- as.matrix(read.delim(ex("Data1.txt"), sep = " ", header = FALSE))
data2 <- as.matrix(read.delim(ex("Data2.txt"), sep = " ", header = FALSE))
mu    <- as.numeric(read.delim(ex("Mean0.txt"), sep = " ", header = FALSE)[, 1])

## Pure-R reference
out_R   <- IWT2(data1, data2, mu, B = 1000)

## Compiled RcppEigen + OpenMP engine
out_cpp <- fastIWT2(data1, data2, mu, B = 1000, MPthreads = 4)

str(out_cpp)
