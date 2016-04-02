library(MASS)
library(VGAM)
library(FactoMineR)

#####
## setup
#####
numObs = 100
numVar = 1000

## create data k = 2 
u = mvrnorm(numObs, mu = rep(0, 2), Sigma = matrix(c(1, 0, 0, 1), nrow = 2) )
tmp = qr(u)
orthU = qr.Q(tmp, complete=FALSE)

v = mvrnorm(numVar, mu = rep(0, 2), Sigma = matrix(c(1, 0, 0, 1), nrow = 2) )
tmp = qr(v)
orthV = qr.Q(tmp, complete=FALSE)

xD = diag(c(4, 3), 2)

#xSig = orthX %*% xD %*% t(orthX)
y = orthU %*% xD %*% t(orthV)
y = y + mvrnorm(numObs, mu = rep(0, numVar), Sigma = diag(mean(diag(cov(y)))/2, numVar))
#test = svd(y)
#y = mvrnorm(numObs, mu = rep(0, 1000), Sigma = xSig) + 
#    mvrnorm(numObs, mu = rep(0, 1000), Sigma = diag(mean(diag(xSig)), numVar))

rSqs = rep(NA, numVar)
for(a in 1:numVar){
    rSqs[a] = summary(lm(y[,a] ~ orthU))$adj.r.squared
}
summary(rSqs)

#fullFac = PCA(y, graph = F)
#fullFac2 = dimdesc(fullFac, axes = 1:2, proba = 1)
#hist(abs(fullFac2$Dim.1$quanti[,1]))
#hist(abs(fullFac2$Dim.2$quanti[,1]))

#foobar = cbind(abs(fullFac2$Dim.1$quanti[order(rownames(fullFac2$Dim.1$quanti)), 1]),
#           abs(fullFac2$Dim.2$quanti[order(rownames(fullFac2$Dim.2$quanti)), 1]))
#foobar = foobar[order(rowMeans(foobar), decreasing = T), ]
#View(foobar)

#orderY = as.numeric(sub("V", "", row.names(foobar)))
#orderY = as.numeric(sub("V", "", row.names(fullFac2$Dim.1$quanti)))
#newY = y[, orderY]
newY = y[, order(rSqs, decreasing = T)]

#####
## newest sims
#####
k = 2
p = c(5, 25, 50, 101, 500, 1000)
epsilon = seq(0.1, 1, 0.1)
numSim = 1
totOut = length(p) * length(epsilon) * numSim
outList = vector("list", length(p))


for(c in 1:length(p)){
    time = proc.time()
    outList[[c]] = vector("list", length(epsilon))
    
    DF = scale(newY[, 1:p[c]])
    #DF = scale(newY[, sample(1:1000, p[c])])
    d = min(p[c], numObs)
    xsvd = svd(DF)
    resid = DF - xsvd$u[, 1:k] %*% diag(xsvd$d[1:k], k) %*% t(xsvd$v[, 1:k])
    
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
                newD[a] = xsvd$d[a] + rlaplace(1, 0, gsL / epsilon[b])
            }
            
            c1 = solve(t(newU) %*% newU)
            c2 =  eigen(c1)$vectors %*% diag(sqrt(eigen(c1)$values), k) %*% t(eigen(c1)$vectors)
            newU = newU %*% c2
            
            DFE = scale(newU %*% diag(newD, k) %*% t(xsvd$v[, 1:2])) + resid[sample(1:numObs, numObs),]
            DFE = scale(newU %*% diag(newD, k) %*% t(xsvd$v[, 1:2]))
            
            boo = PCA(DF, graph = T)
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

resultsMu1 = matrix(NA, ncol = 6, nrow = 10)
resultsMu2 = matrix(NA, ncol = 6, nrow = 10)
for(c in 1:10){
    for(b in 1:6){
        hold = outList[[b]][[c]]
        resultsMu1[c, b] = mean(abs(hold[, 1] - hold[, 3]))
        resultsMu2[c, b] = mean(abs(hold[, 2] - hold[, 4]))
    }
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
## interpret results
#####

## one set
dim(test)
resultsMu1 = matrix(NA, ncol = 6, nrow = 10)
resultsMu2 = matrix(NA, ncol = 6, nrow = 10)
resultsSd1 = matrix(NA, ncol = 6, nrow = 10)
resultsSd2 = matrix(NA, ncol = 6, nrow = 10)
for(c in 1:10){
    for(b in 1:6){
        resultsMu1[c, b] = 0
        resultsMu2[c, b] = 0
        resultsSd1[c, b] = 0
        resultsSd2[c, b] = 0
        for(a in 1:100){
                hold = test[c, (a + (b - 1) * 100)][[1]]
                resultsMu1[c, b] = (mean(abs(hold[, 1] - hold[, 3])) + resultsMu1[c, b] * (a - 1)) / a
                
                resultsMu2[c, b] = (mean(abs(hold[, 2] - hold[, 4])) + resultsMu2[c, b] * (a - 1)) / a
                
                resultsSd1[c, b] = (sd(abs(hold[, 1] - hold[, 3])) + resultsSd1[c, b] * (a - 1)) / a
                
                resultsSd2[c, b] = (sd(abs(hold[, 2] - hold[, 4])) + resultsSd2[c, b] * (a - 1)) / a
        }
    }
}

par(mfrow = c(1, 2))
plot(resultsMu1[,1], ylim = range(resultsMu1), type = "b", col = 1)
for(a in 2:6){
    lines(resultsMu1[, a], type = "b", col = a)
}

plot(resultsMu2[,1], ylim = range(resultsMu2), type = "b", col = 1)
for(a in 2:6){
    lines(resultsMu2[, a], type = "b", col = a)
}

plot(resultsSd1[,1], ylim = range(resultsSd1), type = "b", col = 1)
for(a in 2:6){
    lines(resultsSd1[, a], type = "b", col = a)
}

plot(resultsSd2[,1], ylim = range(resultsSd2), type = "b", col = 1)
for(a in 2:6){
    lines(resultsSd2[, a], type = "b", col = a)
}


## no bootstrap
dim(test2)
results2Mu1 = matrix(NA, ncol = 6, nrow = 10)
results2Mu2 = matrix(NA, ncol = 6, nrow = 10)
results2Sd1 = matrix(NA, ncol = 6, nrow = 10)
results2Sd2 = matrix(NA, ncol = 6, nrow = 10)
for(c in 1:10){
    for(b in 1:6){
        results2Mu1[c, b] = 0
        results2Mu2[c, b] = 0
        results2Sd1[c, b] = 0
        results2Sd2[c, b] = 0
        for(a in 1:100){
            hold = test2[c, (a + (b - 1) * 100)][[1]]
            results2Mu1[c, b] = (mean(abs(hold[, 1] - hold[, 3])) + results2Mu1[c, b] * (a - 1)) / a
            
            results2Mu2[c, b] = (mean(abs(hold[, 2] - hold[, 4])) + results2Mu2[c, b] * (a - 1)) / a
            
            results2Sd1[c, b] = (sd(abs(hold[, 1] - hold[, 3])) + results2Sd1[c, b] * (a - 1)) / a
            
            results2Sd2[c, b] = (sd(abs(hold[, 2] - hold[, 4])) + results2Sd2[c, b] * (a - 1)) / a
        }
    }
}

par(mfrow = c(1, 2))
plot(results2Mu1[,1], ylim = range(results2Mu1), type = "b", col = 1)
for(a in 2:6){
    lines(results2Mu1[, a], type = "b", col = a)
}

plot(results2Mu2[,1], ylim = range(results2Mu2), type = "b", col = 1)
for(a in 2:6){
    lines(results2Mu2[, a], type = "b", col = a)
}

plot(results2Sd1[,1], ylim = range(results2Sd1), type = "b", col = 1)
for(a in 2:6){
    lines(results2Sd1[, a], type = "b", col = a)
}

plot(results2Sd2[,1], ylim = range(results2Sd2), type = "b", col = 1)
for(a in 2:6){
    lines(results2Sd2[, a], type = "b", col = a)
}


## multiple imputes
dim(test3)
results3Mu1 = matrix(NA, ncol = 6, nrow = 10)
results3Mu2 = matrix(NA, ncol = 6, nrow = 10)
results3Sd1 = matrix(NA, ncol = 6, nrow = 10)
results3Sd2 = matrix(NA, ncol = 6, nrow = 10)
for(c in 1:10){
    for(b in 1:6){
        results3Mu1[c, b] = 0
        results3Mu2[c, b] = 0
        results3Sd1[c, b] = 0
        results3Sd2[c, b] = 0
        for(a in 1:100){
            hold = test3[c, (a + (b - 1) * 100)][[1]]
            
            tempMu1 = tempMu2 = tempSd1 = tempSd2 = 0
            for(m in 1:5){
                tempMu1 = (mean(abs(hold[, 1] - hold[, (1 + 2 * m)])) + tempMu1 * (m - 1)) / m
                tempMu2 = (mean(abs(hold[, 2] - hold[, (2 + 2 * m)])) + tempMu2 * (m - 1)) / m
                tempSd1 = (sd(abs(hold[, 1] - hold[, (1 + 2 * m)])) + tempSd1 * (m - 1)) / m
                tempSd2 = (sd(abs(hold[, 2] - hold[, (2 + 2 * m)])) + tempSd2 * (m - 1)) / m
            }
            
            results3Mu1[c, b] = (tempMu1 + results3Mu1[c, b] * (a - 1)) / a
            
            results3Mu2[c, b] = (tempMu2 + results3Mu2[c, b] * (a - 1)) / a
            
            results3Sd1[c, b] = (tempSd1 + results3Sd1[c, b] * (a - 1)) / a
            
            results3Sd2[c, b] = (tempSd2 + results3Sd2[c, b] * (a - 1)) / a
        }
    }
}

par(mfrow = c(2, 2))
plot(results3Mu1[,1], ylim = range(results3Mu1), type = "b")
for(a in 2:6){
    lines(results3Mu1[, a], type = "b")
}

plot(results3Mu2[,1], ylim = range(results3Mu2), type = "b")
for(a in 2:6){
    lines(results3Mu2[, a], type = "b")
}

plot(results3Sd1[,1], ylim = range(results3Sd1), type = "b")
for(a in 2:6){
    lines(results3Sd1[, a], type = "b")
}

plot(results3Sd2[,1], ylim = range(results3Sd2), type = "b")
for(a in 2:6){
    lines(results3Sd2[, a], type = "b")
}

#####
## all plot together
#####
par(mfrow = c(3, 2))

plot(resultsMu1[,1], ylim = range(c(resultsMu1, results2Mu1, results3Mu1)), type = "b", col = 1)
for(a in 4:6){
    lines(resultsMu1[, a], type = "b", col = a)
}

plot(resultsMu2[,1], ylim = range(c(resultsMu2, results2Mu2, results3Mu2)), type = "b", col = 1)
for(a in 4:6){
    lines(resultsMu2[, a], type = "b", col = a)
}

plot(results2Mu1[,1], ylim = range(c(resultsMu1, results2Mu1, results3Mu1)), type = "b", col = 1)
for(a in 4:6){
    lines(results2Mu1[, a], type = "b", col = a)
}

plot(results2Mu2[,1], ylim = range(c(resultsMu2, results2Mu2, results3Mu2)), type = "b", col = 1)
for(a in 4:6){
    lines(results2Mu2[, a], type = "b", col = a)
}

plot(results3Mu1[,1], ylim = range(c(resultsMu1, results2Mu1, results3Mu1)), type = "b")
for(a in 4:6){
    lines(results3Mu1[, a], type = "b", col = a)
}

plot(results3Mu2[,1], ylim = range(c(resultsMu2, results2Mu2, results3Mu2)), type = "b")
for(a in 4:6){
    lines(results3Mu2[, a], type = "b", col = a)
}

par(mfrow = c(3, 2))

plot(resultsSd1[,1], ylim = range(c(resultsSd1, results2Sd1, results3Sd1)), type = "b", col = 1)
for(a in 4:6){
    lines(resultsSd1[, a], type = "b", col = a)
}

plot(resultsSd2[,1], ylim = range(c(resultsSd2, results2Sd2, results3Sd2)), type = "b", col = 1)
for(a in 4:6){
    lines(resultsSd2[, a], type = "b", col = a)
}

plot(results2Sd1[,1], ylim = range(c(resultsSd1, results2Sd1, results3Sd1)), type = "b", col = 1)
for(a in 4:6){
    lines(results2Sd1[, a], type = "b", col = a)
}

plot(results2Sd2[,1], ylim = range(c(resultsSd2, results2Sd2, results3Sd2)), type = "b", col = 1)
for(a in 4:6){
    lines(results2Sd2[, a], type = "b", col = a)
}

plot(results3Sd1[,1], ylim = range(c(resultsSd1, results2Sd1, results3Sd1)), type = "b")
for(a in 4:6){
    lines(results3Sd1[, a], type = "b", col = a)
}

plot(results3Sd2[,1], ylim = range(c(resultsSd2, results2Sd2, results3Sd2)), type = "b")
for(a in 4:6){
    lines(results3Sd2[, a], type = "b", col = a)
}



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





