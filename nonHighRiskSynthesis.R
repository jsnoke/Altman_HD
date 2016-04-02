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

p = c(5, 10, 25, 50, 75, 90)
cc = 1
k = 2

#####
## synthesis
#####

synCompare = function(cc, k, p, newY, numObs){
    DF = scale(newY[, 1:p[cc]])
    xsvd = svd(DF)
    resid = DF - xsvd$u[, 1:k] %*% diag(xsvd$d[1:k], k) %*% t(xsvd$v[, 1:k])
    
    ## normal
    tradSyn = syn(DF, m = 1, method = "cart")
    
    ## svd
    DFE = xsvd$u[, 1:k] %*% diag(xsvd$d[1:k], k) %*% t(xsvd$v[, 1:k]) + resid[sample(1:numObs, numObs),]
    
    ## compare
    trueFac = PCA(DF, graph = F)
    tradFac = PCA(tradSyn$syn, graph = F)
    svdFac = PCA(DFE, graph = F)
    
    outList = cbind(
        dimdesc(trueFac, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
            order(rownames(dimdesc(trueFac, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
        dimdesc(trueFac, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
            order(rownames(dimdesc(trueFac, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),],
        dimdesc(tradFac, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
            order(rownames(dimdesc(tradFac, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
        dimdesc(tradFac, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
            order(rownames(dimdesc(tradFac, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),],
        dimdesc(svdFac, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
            order(rownames(dimdesc(svdFac, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
        dimdesc(svdFac, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
            order(rownames(dimdesc(svdFac, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),])
    
    return(outList)

}

test = synCompare(cc, k, p, newY, numObs)
test


