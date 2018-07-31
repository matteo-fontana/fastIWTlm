IWT2 = function (data1, data2, mu = 0, B = 1000, maxrow = 1, paired = F,
                 recycle = F, alternative = 'two.sided') {
    
    # Alternative
    possible_alternatives = c('two.sided', 'greater', 'less')
    if (!(alternative %in% possible_alternatives))
        stop('Invalid alternative. Choose among ',
             paste0(possible_alternatives, collapse = ', '), '.')
    
    # Data
    if (is.matrix(data1) & is.matrix(data2)) {
        coeff1 = data1
        coeff2 = data2
    } else
        stop ('data1 and data2 must be matrices.')
    
    # Dimensions
    n1 = dim(coeff1)[1]
    n2 = dim(coeff2)[1]
    p  = dim(coeff1)[2]
    n  = n1 + n2
    coeff1 = coeff1 - matrix(data = mu, nrow = n1, ncol = p, byrow = T)
    data.eval = coeff = rbind(coeff1, coeff2)
    data.eval[1:n1,] = data.eval[1:n1,] + matrix(data = mu, nrow = n1, ncol = p, byrow = T)
    
    # print('Point-wise tests')
    
    # Alternative choice
    meandiff = colMeans(coeff[1:n1,]) - colMeans(coeff[(n1+1):n,])
    sign.diff = sign(meandiff)
    sign.diff[which(sign.diff==-1)] = 0
    
    # Test statistic on original data
    T0 = switch(alternative,
                two.sided = (meandiff)^2,
                greater   = (meandiff* sign.diff   )^2,
                less      = (meandiff*(sign.diff-1))^2)
    
    # Test statistic on permuted data
    T_coeff = matrix(ncol = p, nrow = B)
    for (perm in 1 : B) {
        if (paired) {
            if.perm = rbinom(n1, 1, 0.5) 
            coeff_perm = coeff
            for (couple in 1 : n1)
                if (if.perm[couple]==1)
                    coeff_perm[c(couple,n1+couple),] = coeff[c(n1+couple,couple),]
        } else {
            permutazioni = sample(n)
            coeff_perm = coeff[permutazioni,]
        }
        meandiff = colMeans(coeff_perm[1:n1,]) - colMeans(coeff_perm[(n1+1):n,])
        sign.diff = sign(meandiff)
        sign.diff[which(sign.diff==-1)] = 0
        T_coeff[perm,] = switch(alternative,
                                two.sided = (meandiff)^2,
                                greater   = (meandiff* sign.diff   )^2,
                                less      = (meandiff*(sign.diff-1))^2)
    }
    pval = numeric(p)
    for (i in 1:p) pval[i] = sum(T_coeff[,i]>=T0[i]) / B
    
    # print('Interval-wise tests')
    
    # Asymmetric combination matrix
    matrice_pval_asymm = matrix(nrow = p, ncol = p)
    matrice_pval_asymm[p,] = pval[1:p]
    T0_2x = c(T0,T0)
    T_coeff_2x = cbind(T_coeff,T_coeff)
    
    if (recycle==T) {
        for (i in (p-1) : maxrow) { # rows
            for (j in 1 : p) { # columns
                inf = j
                sup = (p-i)+j
                T0_temp = sum(T0_2x[inf:sup])
                T_temp = rowSums(T_coeff_2x[,inf:sup])
                pval_temp = sum(T_temp>=T0_temp) / B
                matrice_pval_asymm[i,j] = pval_temp
            }
            #print(paste('creating the p-value matrix: end of row ',
            #      as.character(p-i+1), ' out of ', as.character(p), sep = ''))
        }
    } else
        for (i in (p-1) : maxrow) { # rows
            for (j in 1 : i) { # columns
                inf = j
                sup = (p-i)+j
                T0_temp = sum(T0_2x[inf:sup])
                T_temp = rowSums(T_coeff_2x[,inf:sup])
                pval_temp = sum(T_temp>=T0_temp) / B
                matrice_pval_asymm[i,j] = pval_temp
            }
            #print(paste('creating the p-value matrix: end of row ',
            #      as.character(p-i+1), ' out of ', as.character(p), sep = ''))
        }
    
    pval.correct = function (pval.matrix) {
        matrice_pval_2_2x = cbind(pval.matrix, pval.matrix)
        p = dim(pval.matrix)[2]
        matrice_pval_2_2x = matrice_pval_2_2x[,(2*p):1]
        corrected.pval = numeric(p)
        corrected.pval.matrix = matrix(nrow = p, ncol = p)
        corrected.pval.matrix[p,] = pval.matrix[p,p:1]
        
        for (var in 1 : p) {
            pval_var = matrice_pval_2_2x[p,var]
            inizio = var
            fine = var
            for (riga in (p-1) : 1) {
                fine = fine + 1
                pval_cono = matrice_pval_2_2x[riga, inizio:fine]
                pval_var = max(pval_var, pval_cono,na.rm = T)
                corrected.pval.matrix[riga,var] = pval_var
            }
            corrected.pval[var] = pval_var
        }
        
        corrected.pval = corrected.pval[p:1]
        corrected.pval.matrix = corrected.pval.matrix[maxrow,p:1]
        return(corrected.pval.matrix)
    }
    
    #print('Interval-wise testing completed!')
    
    corrected.pval.matrix = pval.correct(matrice_pval_asymm)
    
    return (list (T0 = T0,
                  point = pval,
                  inter = matrice_pval_asymm,
                  corre = corrected.pval.matrix) )
}
