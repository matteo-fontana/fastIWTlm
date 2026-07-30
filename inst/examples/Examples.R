## Load auxiliary functions
source('AuxiliaryFunctions.R')

## Define input parameters
n1 = n2 = 30
p = 200
B = 2000
alt = 'two.sided'
maxrow = 0
paired = T
recycle = F
THREADS = 1

## Generate random data
library(mvtnorm)
data1 = ComputeData(n = n1, p = p, mu1 = 0, mu2 = 2)
data2 = ComputeData(n = n2, p = p, mu1 = 0, mu2 = 0)
mu = rep(0,p)

## Plot input data
PlotData(data1, data2)

## Write text files
write.table(c(n1,n2,p), file = 'Param.txt', col.names = F, row.names = F)
write.table(data1, file = 'Data1.txt', col.names = F, row.names = F)
write.table(data2, file = 'Data2.txt', col.names = F, row.names = F)
write.table(mu,    file = 'Mean0.txt', col.names = F, row.names = F)

## Load R and C++ versions
source('IWT2.R')
library(Rcpp)
Sys.setenv('PKG_CXXFLAGS' = '-std=c++11 -Wall -pedantic -fopenmp')
sourceCpp('IWT2OMP.cpp')

## First example
maxrow = 0
recycle = F
RES = Execute(n1, n2, p, B, alt, maxrow, paired, recycle, THREADS) # 50 seconds
PlotPval(RES)
PlotImage(RES, data1, data2, mu)

## Second example
maxrow = p/2
recycle = T
RES = Execute(n1, n2, p, B, alt, maxrow, paired, recycle, THREADS) # 30 seconds
PlotPval(RES)
PlotImage(RES, data1, data2, mu)
