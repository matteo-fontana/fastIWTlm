GeneraDati = function(n1 = 100, n2 = 100, p = 200) {
    
    # Loading libraries
    library(mvtnorm)
    source('AuxiliaryFunctions.R')
    
    # Generating data
    data1 = ComputeData(n = n1, p = p, mu1 = 0, mu2 = 2)
    data2 = ComputeData(n = n2, p = p, mu1 = 0, mu2 = 0)
    mu = c(rep(0,p/2),rep(0,p/2))
    
    # Writing files
    write.table(c(n1,n2,p), file = 'Param.txt', col.names = F, row.names = F)
    write.table(data1, file = 'Data1.txt', col.names = F, row.names = F)
    write.table(data2, file = 'Data2.txt', col.names = F, row.names = F)
    write.table(mu,    file = 'Mean0.txt', col.names = F, row.names = F)
}
