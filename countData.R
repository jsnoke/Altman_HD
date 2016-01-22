#######################

##Count data

#######################
library(matrixStats)
library(MASS)
library(synthpop)
library(FactoMineR)

nuU = rnorm(numObs, mean = 1, sd = 1)
nuZ = mvrnorm(100, mu = rep(0, 5), Sigma = diag(0.1, 5))

nuX = nuU + nuZ
lambda = exp(nuX)    
summary(lambda)

xDat = data.frame(matrix(NA, nrow = 100, ncol = 5))
for(i in 1:100){
    for(j in 1:5){
        xDat[i,j] = rpois(1, lambda[i,j])   
    }
}
summary(xDat)

test = PCA(xDat, scale.unit = F)
testNu = PCA(nuX, scale.unit = F)
newNu = reconst(testNu, ncp = 1)
newLambda = exp(newNu)

newXDat = data.frame(matrix(NA, nrow = 100, ncol = 5))
for(i in 1:100){
    for(j in 1:5){
        newXDat[i,j] = rpois(1, newLambda[i,j])   
    }
}

synX = syn(xDat, m = 5)
summary(synX$syn[[1]])
View(synX$syn[[1]])

#---------------------------
##two factors

nuU1 = mvrnorm(numObs, mu = rep(0.5, 2), Sigma = diag(0.5, 2))
nuZ1 = mvrnorm(100, mu = rep(0, 5), Sigma = diag(0.1, 5))

nuX1 = cbind(nuU1[, 1] + nuZ1[, 1:2], nuU1[, 2] + nuZ1[, 3:5])
lambda1 = exp(nuX1)    
summary(lambda1)
xDat1 = data.frame(matrix(NA, nrow = 100, ncol = 5))
for(i in 1:100){
    for(j in 1:5){
        xDat1[i,j] = rpois(1, lambda1[i,j])   
    }
}


##synthetic by PC actual
PC = PCA(nuX1, scale.unit = F)
newNu1 = reconst(PC, ncp = 2)

residDat = nuX1 - newNu1

synSamp = matrix(NA, ncol = 5, nrow = numObs)
bootResid = vector("list", 5)
secondSamp = matrix(NA, ncol = 5, nrow = numObs)
secondResid = vector("list", 5)
synDat = vector("list", 5)

for(i in 1:5){
    synSamp[, i] = sample(1:numObs, 100, replace = T)
    secondSamp[, i] = sample(1:numObs, 100, replace = T)
    
    bootResid[[i]] = matrix(NA, ncol = 5, nrow = numObs)
    secondResid[[i]] = matrix(NA, ncol = 5, nrow = numObs)
    for(j in 1:numObs){
        bootResid[[i]][j, ] = residDat[synSamp[j,i], ] 
        secondResid[[i]][j, ] = residDat[secondSamp[j,i], ] 
    }
    
    temp = newNu1 + bootResid[[i]]
    tempPC = PCA(temp, scale.unit = F)
    tempRecon = reconst(tempPC, ncp = 2)
    
    synDat[[i]] = tempRecon + secondResid[[i]]
    colnames(synDat[[i]]) = paste("x", 1:5 , sep = "")
    
}

newLambda1 = vector("list", 5)
for(i in 1:5){
    newLambda1[[i]] = exp(synDat[[i]])
}

newXDat1 = vector("list", 5)
for(k in 1:5){
    newXDat1[[k]] =  data.frame(matrix(NA, nrow = 100, ncol = 5))
    for(i in 1:100){
        for(j in 1:5){
            newXDat1[[k]][i,j] = rpois(1, newLambda1[[k]][i,j])   
        }
    }
}

synX2 = syn(xDat1, m = 5, method = "parametric")
summary(synX$syn[[2]])














