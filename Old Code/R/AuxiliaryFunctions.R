# Run both versions of the algorithm with the same inputs
Execute = function (n1 = 50, n2 = 50, p = 100, B = 1000, alt = 'two.sided',
                    maxrow = 0, paired = T, recycle = F, THREADS = 1) {
    
    ## Print information about the current execution
    cat('\n****************************************************************\n')
    cat(paste('n','p','B','maxr','pair','recy','alternative','THRE', sep='\t'))
    cat('\n')
    cat(paste(n1, p, B, maxrow, paired, recycle, alt, THREADS, sep = '\t'))
    cat('\n')
    
    ## Execute old version of the algorithm
    T1 = Sys.time()
    IWT2(as.matrix(read.table('Data1.txt')),
         as.matrix(read.table('Data2.txt')),
         as.matrix(read.table('Mean0.txt')),
         B, maxrow+1, paired, recycle, alt)
    T2 = Sys.time()
    cat('\nOld version (R): ')
    print(T2-T1)
    
    ## Execute new version of the algorithm
    T2 = Sys.time()
    RES = IWT2OMP(B, alt, maxrow, paired, recycle, THREADS)
    T3 = Sys.time()
    cat('\nNew version (C++): ')
    print(T3-T2)
    return(RES)
}

# Generate random matrix with different mean values in the two halves
ComputeData = function (n = 50, p = 100, mu1 = 0, mu2 = 0) {
    mu = c(rep(mu1, floor(p/2)), rep(mu2, ceiling(p/2)))
    return(rmvnorm(n, mean = mu, sigma = diag(p)/4))
}

# Plot input data and their mean difference
PlotData = function (data1, data2) {
    x11()
    par(mfrow = c(2,1))
    matplot(t(data1), type = 'l', col = 'blue', main = 'Input data')
    matplot(t(data2), type = 'l', col = 'red', add = T)
    legend('topright', c('data1','data2'), col = c('blue','red'), lty = c(2,2))
    matplot(colMeans(data1) - colMeans(data2), main = 'Mean difference',
            type = 'l', col = 'black')
}

# Plot unadjusted and adjusted p-values
PlotPval = function (RES) {
    x11()
    par(mfrow = c(2,1))
    plot(RES$point, main = 'Unadjusted p-values', ylab = 'p-value')
    plot(RES$corre, main = 'Adjusted p-values', ylab = 'p-value')
}

# Plot heatmap, adjusted p-values and input data
PlotImage = function (RES, data1, data2, mu) {
    source('IWTimage.R')
    result = list(test = '2pop',
                  mu = mu,
                  corre = RES$corre,
                  point = RES$point,
                  inter = RES$inter,
                  data.eval = rbind(data1 - matrix(mu, nrow(data1), ncol(data1), byrow = T),
                                    data2),
                  ord_labels = c(rep(1, nrow(data1)), rep(2, nrow(data2))))
    class(result) = 'IWT2'
    x11()
    IWTimage(result)
}
