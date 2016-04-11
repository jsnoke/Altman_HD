library(MASS)
library(synthpop)
library(FactoMineR)

numObs = 100
numVar = 1000

## create data k = 2 
u = mvrnorm(numObs, mu = rep(0, 2), Sigma = matrix(c(1, 0, 0, 1), nrow = 2) )
tmp = qr(u)
orthU = qr.Q(tmp, complete=FALSE)

v = mvrnorm(numVar, mu = rep(0, 2), Sigma = matrix(c(1, 0, 0, 1), nrow = 2) )
tmp = qr(v)
orthV = qr.Q(tmp, complete=FALSE)

xD = diag(c(4, 2), 2)

y = orthU %*% xD %*% t(orthV)
y = y + mvrnorm(numObs, mu = rep(0, numVar), Sigma = diag(mean(diag(cov(y)))/2, numVar))

rSqs = rep(NA, numVar)
for(a in 1:numVar){
    rSqs[a] = summary(lm(y[,a] ~ orthU))$r.squared
}
summary(rSqs)

newY = y[, order(rSqs, decreasing = T)]

p = c(5, 25, 50, 101, 500, 1000)
cc = 1
k = 2
numSim = 100

#########
## new risk eval
#########

rowDiff = function(df, row){
    return(abs(row - df))
}

closest = function(list, tol){
    return(which(rowSums((list < tol)) == max(rowSums((list < tol)))))
}


matchesT = vector("list", length(p))
matchMSET = vector("list", length(p))
matchesS = vector("list", length(p))
matchMSES = vector("list", length(p))

####
for(cc in 1:length(p)){
    matchesT[[cc]] = vector("list", numSim)
    matchMSET[[cc]] = rep(NA, numSim)
    matchesS[[cc]] = vector("list", numSim)
    matchesS[[cc]] = vector("list", numSim)
    matchMSES[[cc]] = rep(NA, numSim)
    
    DF = scale(newY[, 1:p[cc]])
    xsvd = svd(DF)
    resid = DF - xsvd$u[, 1:k] %*% diag(xsvd$d[1:k], k) %*% t(xsvd$v[, 1:k])
    
    for(n in 1:numSim){
        matchesS[[cc]][[n]] = matrix(NA, nrow = numObs, ncol = 2)
        ####
        ### normal
        tradSyn = syn(DF, m = 1, method = "cart")$syn
        
        foo = t(apply(as.matrix(tradSyn[, 1:(p[cc] - 1)]), 1, rowDiff, 
                      row = DF[1, 1:(p[cc] - 1), drop = F]))
        
        closeInd = closest(foo, sd(DF))
        matchesT[[cc]][[n]] = as.numeric(rowSums(foo < sd(DF))[closeInd])[1]
        
        allMatch =  abs(tradSyn[closeInd, p[cc]] - DF[1, p[cc]])
        
        if(is.na(var(allMatch)) == TRUE){
            matchMSET[[cc]][n] = mean(allMatch)^2
        } else
            matchMSET[[cc]][n] = mean(allMatch)^2 + var(allMatch)
        
        ####
        ## svd
        DFE = xsvd$u[, 1:k] %*% diag(xsvd$d[1:k], k) %*% t(xsvd$v[, 1:k]) + resid[sample(1:numObs, numObs),]
        
        matchesS[[cc]][[n]] = numMatch(DFE, DF, epsilon = 0.1)
        foo = t(apply(as.matrix(DFE[, 1:(p[cc] - 1)]), 1, rowDiff, 
                      row = DF[1, 1:(p[cc] - 1), drop = F]))
        
        closeInd = closest(foo, sd(DF))
        matchesS[[cc]][[n]] = as.numeric(rowSums(foo < sd(DF))[closeInd])[1]
        
        allMatch =  abs(DFE[closeInd, p[cc]] - DF[1, p[cc]])
        
        if(is.na(var(allMatch)) == TRUE){
            matchMSES[[cc]][n] = mean(allMatch)^2
        } else
            matchMSES[[cc]][n] = mean(allMatch)^2 + var(allMatch)
    
    }
    cat(cc, '\n')
}



#####
#####
#####
numMatch = function(syn, real, epsilon = 0.1){
    matches = matrix(NA, nrow = nrow(real), ncol = nrow(syn))
    for(i in 1:nrow(syn)){
        for(j in 1:nrow(real)){
            matches[j, i] = sum(abs(syn[i,] - real[j,]) < epsilon)
        }
    }
    
    amount = matrix(NA, nrow = nrow(syn), ncol = 2)
    for(i in 1:ncol(matches)){
        amount[i, 1] = max(matches[, i])
        amount[i, 2] = 1 / sum(matches[, i] == max(matches[, i]))
    }
    return(amount)
}

for(cc in 1:length(p)){

    matchesS[[cc]] = vector("list", numSim)
    
    DF = scale(newY[, 1:p[cc]])
    xsvd = svd(DF)
    resid = DF - xsvd$u[, 1:k] %*% diag(xsvd$d[1:k], k) %*% t(xsvd$v[, 1:k])
    
    for(n in 1:numSim){
        matchesS[[cc]][[n]] = matrix(NA, nrow = numObs, ncol = 2)
        
        ####
        ####
        ## svd
        DFE = xsvd$u[, 1:k] %*% diag(xsvd$d[1:k], k) %*% t(xsvd$v[, 1:k]) + resid[sample(1:numObs, numObs),]
        
        matchesS[[cc]][[n]] = numMatch(DFE, DF, epsilon = 0.1)
        
    }
    
    cat(cc, '\n')
}

totMatch = matrix(NA, nrow = 100, ncol = 6)
summaryMatch = matrix(NA, nrow = 100, ncol = 6)
totMatch01 = matrix(NA, nrow = 100, ncol = 6)
summaryMatch01 = matrix(NA, nrow = 100, ncol = 6)
totMatch05 = matrix(NA, nrow = 100, ncol = 6)
summaryMatch05 = matrix(NA, nrow = 100, ncol = 6)

for(a in 1:6){
    for(b in 1:100){
        totMatch01[b, a] = mean(matchesS[[a]][[b]][ , 1])
        summaryMatch01[b, a] = mean(matchesS[[a]][[b]][ , 2])
    }
}

nonHighRisk = list(totMatch, totMatch01, totMatch05, summaryMatch, summaryMatch01, summaryMatch05)
save(nonHighRisk, file = "nonHighRisk.RData")



