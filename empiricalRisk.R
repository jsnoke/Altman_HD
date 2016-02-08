#######################

##Risk measure - empirical distribution

#######################

library(MASS)
library(synthpop)

numObs = 100
numVar = 5

u = rnorm(numObs, mean = 2, sd = 1)
b = mvrnorm(numObs, mu = rep(0, numVar), Sigma = diag(0.1, numVar))
x = u + b

## function
empirRisk = function(syn, real){
    empiry = matrix(NA, nrow = nrow(real), ncol = nrow(syn))
    for(i in 1:nrow(syn)){
        postProb = matrix(NA, ncol = ncol(real), nrow = nrow(real))
        for(j in 1:ncol(real)){
            diff = abs(syn[i, j] - real[, j])
            temp = ecdf(diff)
            postProb[, j] = (1 + 1/nrow(real)) - temp(diff)
        }
        for(z in 1:nrow(postProb)){
            empiry[z, i] = prod(postProb[z, ]) 
        }
    }
    return(empiry)
}

test = empirRisk(x, x)
apply(test, 2, max)/colSums(test)


##simulation
eps = seq(from = 0.1, to = 5, by = 0.1)
plotfram = matrix(NA, nrow = length(eps), ncol = 6)
for(i in 1:length(eps)){
    
    ## noise
    xsvd = svd(x) 
    tempD = rep(NA, ncol(x))
    tempU = matrix(NA, nrow = nrow(x), ncol = ncol(x))
    tempD[1] = xsvd$d[1]
    tempU[, 1] = xsvd$u[, 1]
    for(z in 2:numVar){
        tempD[z] = xsvd$d[z] + rexp(1, rate = xsvd$d[z]/eps[i])
        tempU[, z] = xsvd$u[, z] + runif(1, min = -eps[i], max = eps[i])
    }
    xE = (tempU %*% diag(tempD) %*% t(xsvd$v))
    
    out = empirRisk(xE, x)
    out2 = apply(out, 2, max)
    plotfram[i, ] = c(median(out2), mean(out), max(out2))
}

par(mfrow = c(1,3))
plot(plotfram[, 1])
plot(plotfram[, 2])
plot(plotfram[, 3])




