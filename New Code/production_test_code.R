library(Rcpp)
Sys.setenv('PKG_CXXFLAGS' = '-std=c++11 -Wall -pedantic -fopenmp')
sourceCpp('New Code/fastIWT2.cpp')

load('New Code/test-data.rdata')

###directly call function

results=fastIWT2(data1,data2,mu,B=5000,MPthreads=4)
