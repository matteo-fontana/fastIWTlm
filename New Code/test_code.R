library(Rcpp)
Sys.setenv('PKG_CXXFLAGS' = '-std=c++11 -Wall -pedantic -fopenmp')
sourceCpp('New Code/IWT2C_exp.cpp')

output=IWT2C_exp()


#now, try input via R
##prepare input objects: first step with known n1-n2-p

data1=as.matrix(read.delim('Data1.txt',sep=' ',header=F))
data2=as.matrix(read.delim('Data2.txt',sep=' ',header=F))
mu=as.matrix(read.delim('Mean0.txt',sep=' ',header=F))

n1=dim(data1)[1]
n2=dim(data2)[1]
p=dim(data2)[2]

save(data1,data2,mu, file='New Code/test-data.rdata')

#first, only input parameters

output=IWT2C_exp(data1,data2,mu)

##experiment, read directly as Eigen matrices






