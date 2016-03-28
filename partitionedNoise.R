library(MASS)
library(VGAM)
library(FactoMineR)

#####
# ignore
#####
x = mvrnorm(100, mu = rep(0, 6), Sigma = diag(1, 6))
temp = qr(x)
orthX = qr.Q(temp, complete = FALSE)

x1 = orthX[, 1:3]
x2 = orthX[, 4:6]

x1E = x1 + x1 %*% mvrnorm(3, mu = rep(0, 3), Sigma = diag(1, 3))

x2E = x2 + x2 %*% mvrnorm(3, mu = rep(0, 3), Sigma = diag(1, 3))

xsvd = svd(x)
hist(x -  xsvd$u[, 1:2, drop = F] %*% diag(xsvd$d, 2) %*% t(xsvd$v[, 1:2]))

rx1 = xsvd$u[, 1, drop = F] %*% diag(xsvd$d, 1) %*% t(xsvd$v[, 1])
rx2 = xsvd$u[, 2:6, drop = F] %*% diag(xsvd$d, 5) %*% t(xsvd$v[, 2:6])

foo = PCA(x, graph = F)

#####
## setup
#####
numObs = 100
numVar = 1000
beta = c(0.5, 3)
epsilon = 0.5

## create data k = 1
x = rnorm(numObs, mean = 2, sd = 1)
y = matrix(NA, nrow = numObs, ncol = numVar)
for(i in 1:ncol(y)){
    y[, i] = beta[1] + beta[2]*x + rnorm(numObs)
}
## create data k = 2 
x = mvrnorm(numObs, mu = rep(2, 2), Sigma = matrix(c(3, 0, 0, 3), nrow = 2) )
y = matrix(NA, nrow = numObs, ncol = numVar)
for(i in 1:ncol(y)){
    if(runif(1) > 0.5)
        y[, i] = beta[1] + beta[2]*x[, 1] + rnorm(numObs)
    else
        y[, i] = beta[2] + beta[1]*x[, 2] + rnorm(numObs)
}

#####
## newest sims
#####
rsqs = matrix(NA, nrow = 1000, ncol = 1)
for(a in 1:1000){
    rsqs[a, 1] = summary(lm(y[ , a] ~ x))$adj.r.squared
}
summary(rsqs)
orderY = y[, order(rsqs, decreasing = T)]

k = 2
p = c(5, 25, 50, 101, 500, 1000)
epsilon = seq(0.1, 1, 0.1)
numSim = 10
totOut = length(p) * length(epsilon) * numSim
outList = vector("list", length(p))


for(c in 1:length(p)){
    time = proc.time()
    outList[[c]] = vector("list", length(epsilon))
    
    DF = scale(orderY[, 1:p[c]])
    d = min(p[c], numObs)
    xsvd = svd(DF)
    resid = DF - xsvd$u[, 1:2] %*% diag(xsvd$d[1:2], 2) %*% t(xsvd$v[, 1:2])
    
    gsV = sqrt(p[c] * numObs) / (xsvd$d[1] - xsvd$d[2])
    gsL = sqrt(p[c] * k)
    
    for(b in 1:length(epsilon)){
        outList[[c]][[b]] = matrix(NA, nrow = p[c] * numSim, ncol = 4)
        colnames(outList[[c]][[b]]) = c("real1", "real2", "noise1", "noise2")
    
        for(n in 1:numSim){
            newU = matrix(NA, nrow = numObs, ncol = k)
            newD = rep(NA, k)
            for(a in 1:k){
                newU[, a] = xsvd$u[, a, drop = F] + rlaplace(numObs, 0, gsV / epsilon[b])
                newD[a] = xsvd$d[a] + rlaplace(1, 0, gsL)
            }
            
            c1 = solve(t(newU) %*% newU)
            c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
            newU = newU %*% c2
            
            DFE = scale(newU %*% diag(newD, k) %*% t(xsvd$v[, 1:2])) + resid[sample(1:numObs, numObs),]
            
            boo = PCA(DF, graph = F)
            foo = PCA(DFE, graph = F)
            
            outList[[c]][[b]][(1 + (n - 1) * p[c]):((n) * p[c]), 1:4] = cbind(
                dimdesc(boo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
                    order(rownames(dimdesc(boo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
                dimdesc(boo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
                    order(rownames(dimdesc(boo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),],
                dimdesc(foo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
                    order(rownames(dimdesc(foo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
                dimdesc(foo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
                    order(rownames(dimdesc(foo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),])
        }
        cat(b)    
    }
    print(proc.time() - time)
    cat(c, "\n")
    
}


summary(abs(outList[[1]][[10]][, 1] - outList[[1]][[10]][, 3]))
hist(outList[[1]][[10]][, 1])
hist(outList[[1]][[10]][, 3])

#DFE = scale(newU %*% diag(newD, k) %*% t(xsvd$v[, 1:2]))
#summary(simOut)
#summary(simOut[,1] - simOut[,3])
#summary(simOut[,2] - simOut[,4])
#hist(simOut[,1] - simOut[,3])
#hist(simOut[,2] - simOut[,4])

#####
## just add noise and re-orth
#####
k = 2
p = c(5, 25, 50, 101, 500, 1000)
epsilon = seq(0.1, 1, 0.1)
numSim = 100
totOut = length(p) * length(epsilon) * numSim
outList = vector("list", length(p))

for(c in 1:length(p)){
    outList[[c]] = vector("list", length(epsilon))
    
    DF = scale(y[, 1:p[c]])
    d = min(p[c], numObs)
    xsvd = svd(DF)
    
    gsV = sqrt(p[c] * numObs) / (xsvd$d[1] - xsvd$d[2])
    gsL = sqrt(p[c] * k)
    
    for(b in 1:length(epsilon)){
        outList[[c]][[b]] = matrix(NA, nrow = p[c] * numSim, ncol = 4)
        colnames(outList[[c]][[b]]) = c("real1", "real2", "noise1", "noise2")
        
        for(n in 1:numSim){
            newU = matrix(NA, nrow = numObs, ncol = k)
            newV = matrix(NA, nrow = p[c], ncol = k)
            newD = rep(NA, k)
            for(a in 1:k){
                newU[, a] = xsvd$u[, a, drop = F] + rlaplace(numObs, 0, gsV / epsilon[b])
                newV[, a] = xsvd$v[, a, drop = F] + rlaplace(p[c], 0, gsV / epsilon[b])
                newD[a] = xsvd$d[a] + rlaplace(1, 0, gsL)
            }
            
            c1 = solve(t(newV) %*% newV)
            c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
            newV = newV %*% c2
            
            c1 = solve(t(newU) %*% newU)
            c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
            newU = newU %*% c2
            
            DFE = newU %*% diag(newD, k) %*% t(newV)
            
            boo = PCA(DF, graph = F)
            foo = PCA(DFE, graph = F)
            
            outList[[c]][[b]][(1 + (n - 1) * p[c]):((n) * p[c]), 1:4] = 
                                        cbind(dimdesc(boo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F],
                                              dimdesc(boo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F],
                                              dimdesc(foo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F],
                                              dimdesc(foo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])
        }
        
        cat(b)
    }
    
    cat(c, "\n")
}

save(outList, file = "addNoiseSim.RData")

#####
## visualize
#####
par(mfrow = c(2,2))
for(b in 1:6){
    for(a in 1:4){
        hist(outList[[b]][[1]][, a])
    }
}

#####
## add noise, re-orth,just to u
#####

DF = scale(y[, sample(1:1000, p[c], replace = F)])
d = min(p[c], numObs)
xsvd = svd(DF)

gsV = sqrt(p[c] * numObs) / (xsvd$d[1] - xsvd$d[2])
gsL = sqrt(p[c] * k)

newU = matrix(NA, nrow = numObs, ncol = k)
#newV = matrix(NA, nrow = p[c], ncol = k)
newD = rep(NA, k)
for(a in 1:k){
    newU[, a] = xsvd$u[, a, drop = F] + rlaplace(numObs, 0, gsV / epsilon[b])
    #newV[, a] = xsvd$v[, a, drop = F] + rlaplace(p[c], 0, gsV / epsilon[b])
    newD[a] = xsvd$d[a] + rlaplace(1, 0, gsL)
}

c1 = solve(t(newV) %*% newV)
c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
newV = newV %*% c2

c1 = solve(t(newU) %*% newU)
c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
newU = newU %*% c2

DFE = scale(newU %*% diag(newD, k) %*% t(xsvd$v[, 1:2]))
#DFE = scale(newU %*% diag(newD, k) %*% t(newV))
#DFE = scale(xsvd$u[, 1:2] %*% diag(newD, k) %*% t(newV))
summary(DFE)

esvd = svd(DFE)

summary(esvd$u)
summary(newU)
summary(xsvd$u)

summary(esvd$v)
summary(newV)
summary(xsvd$v)

esvd$d
newD
xsvd$d

boo = PCA(DF, graph = T)
foo = PCA(DFE, graph = T)


simOut = cbind(dimdesc(boo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
    order(rownames(dimdesc(boo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
    dimdesc(boo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
        order(rownames(dimdesc(boo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),],
    dimdesc(foo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
        order(rownames(dimdesc(foo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
    dimdesc(foo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
        order(rownames(dimdesc(foo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),])
summary(simOut[,1] - simOut[,3])
summary(simOut[,2] - simOut[,4])
hist(simOut[,1] - simOut[,3])
hist(simOut[,2] - simOut[,4])

#####
## add noise, re-orth, synth
#####

DF = scale(y[, sample(1:1000, p[c], replace = F)])
d = min(p[c], numObs)
xsvd = svd(DF)
resid = DF - xsvd$u[, 1:2] %*% diag(xsvd$d[1:2], 2) %*% t(xsvd$v[, 1:2])

gsV = sqrt(p[c] * numObs) / (xsvd$d[1] - xsvd$d[2])
gsL = sqrt(p[c] * k)

newU = matrix(NA, nrow = numObs, ncol = k)
#newV = matrix(NA, nrow = p[c], ncol = k)
newD = rep(NA, k)
for(a in 1:k){
    newU[, a] = xsvd$u[, a, drop = F] + rlaplace(numObs, 0, gsV / epsilon[b])
    #newV[, a] = xsvd$v[, a, drop = F] + rlaplace(p[c], 0, gsV / epsilon[b])
    newD[a] = xsvd$d[a] + rlaplace(1, 0, gsL)
}

#c1 = solve(t(newV) %*% newV)
#c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
#newV = newV %*% c2

c1 = solve(t(newU) %*% newU)
c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
newU = newU %*% c2

DFE = scale(newU %*% diag(newD, k) %*% t(xsvd$v[, 1:2]))
DFE = scale(newU %*% diag(newD, k) %*% t(xsvd$v[, 1:2])) + resid[sample(1:numObs, numObs),]
#DFE = scale(newU %*% diag(newD, k) %*% t(newV))
#DFE = scale(xsvd$u[, 1:2] %*% diag(newD, k) %*% t(newV))
summary(DFE)

esvd = svd(DFE)

boo = PCA(DF, graph = T)
foo = PCA(DFE, graph = T)

simOut = cbind(dimdesc(boo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
    order(rownames(dimdesc(boo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
    dimdesc(boo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
        order(rownames(dimdesc(boo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),],
    dimdesc(foo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F][
        order(rownames(dimdesc(foo, axes = 1:2, proba = 1)$Dim.1$quanti[, 1, drop = F])),],
    dimdesc(foo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F][
        order(rownames(dimdesc(foo, axes = 1:2, proba = 1)$Dim.2$quanti[, 1, drop = F])),])
simOut
View(simOut)
summary(simOut[,1] - simOut[,3])
summary(simOut[,2] - simOut[,4])
hist(simOut[,1] - simOut[,3])
hist(simOut[,2] - simOut[,4])





