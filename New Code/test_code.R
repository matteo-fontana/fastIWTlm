setwd('New Code/')
library(Rcpp)
Sys.setenv('PKG_CXXFLAGS' = '-std=c++11 -Wall -pedantic -fopenmp')
sourceCpp('IWT2OMP.cpp')
