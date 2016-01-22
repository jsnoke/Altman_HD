library(MASS)
library(FactoMineR)

#####
##### Setup
#####
numObs = 100
numVar = 5

u = rnorm(numObs, mean = 2, sd = 1)
b = mvrnorm(numObs, mu = rep(0, numVar), Sigma = diag(0.1, numVar))
x = u + b

xsvd = svd(x)

#####
##### Add noise to d
#####

epsilon = 2

## add noise
temp = rep(NA, 5)
temp[1] = xsvd$d[1]
for(i in 2:numVar){
    temp[i] = xsvd$d[i] + rexp(1, rate = epsilon/xsvd$d[i])
}
xE = (xsvd$u %*% diag(temp) %*% t(xsvd$v))
xEsvd = svd(xE)

## reverse noise
temp2 = rep(NA, 5)
temp2[1] = xEsvd$d[1]
for(i in 2:numVar){
    temp2[i] = xEsvd$d[i] - rexp(1, rate = epsilon/xEsvd$d[i])
}
xRev = (xEsvd$u %*% diag(temp2) %*% t(xEsvd$v))
xRsvd = svd(xRev)

## compare
xsvd$d
xEsvd$d
xRsvd$d

## u and v should be same within sign change
head(xsvd$u)
head(xEsvd$u)
head(xRsvd$u)
summary(xsvd$u[,5] + xEsvd$u[,4])

xsvd$v
xEsvd$v
xRsvd$v

## pca
xPC = PCA(x)
xRevPC = PCA(xRev)
xEPC = PCA(xE)

xPC$eig
xEPC$eig
xRevPC$eig

xPC$var$coord
xEPC$var$coord
xRevPC$var$coord

## risk
summary(x)
summary(xE)
summary(xRev)
summary(abs(x - xE))
summary(abs(x - xRev))



#####
##### Add noise to u
#####

epsilon = 0.25

## add noise
temp = matrix(NA, nrow = numObs, ncol = numVar)
temp[, 1] = xsvd$u[, 1]
for(i in 2:numVar){
    temp[, i] = xsvd$u[, i] + runif(numObs, min = -epsilon, max = epsilon)
}
xE = (temp %*% diag(xsvd$d) %*% t(xsvd$v))
xEsvd = svd(xE)

## reverse noise
temp2 = matrix(NA, nrow = numObs, ncol = numVar)
temp2[, 1] = xEsvd$u[, 1]
for(i in 2:numVar){
    temp2[, i] = xEsvd$u[, i] - runif(numObs, min = -epsilon, max = epsilon)
}
xRev = (temp %*% diag(xEsvd$d) %*% t(xEsvd$v))
xRsvd = svd(xRev)

## compare
xsvd$d
xEsvd$d
xRsvd$d

## 
head(xsvd$u)
head(xEsvd$u)
head(xRsvd$u)
summary(xsvd$u)
summary(xEsvd$u)
summary(xRsvd$u)
summary(xsvd$u[,5] + xEsvd$u[,4])

xsvd$v
xEsvd$v
xRsvd$v

## pca
xPC = PCA(x)
xRevPC = PCA(xRev)
xEPC = PCA(xE)

xPC$eig
xEPC$eig
xRevPC$eig

xPC$var$coord
xEPC$var$coord
xRevPC$var$coord

## utility
summary(abs(xPC$var$coord[,1] - xEPC$var$coord[,1]))

## risk
summary(rowMeans(abs(x - xE)))


summary(x)
summary(xE)
summary(xRev)
summary(abs(x - xE))
summary(abs(x - xRev))


#####
##### Simulation
#####

outRisk = data.frame(matrix(NA, nrow = 300, ncol = 10))
outUtil = data.frame(matrix(NA, nrow = 300, ncol = 10))
#epsilon = c(0.001, 0.025, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 0.75, 1)
#alpha = c(25, 10, 5, 4, 3, 2, 1, 0.5, 0.25, 0.1)
epsilon = c(1.5, 3, 5, 10, 25)
alpha = c(0.05, 0.025, 0.01, 0.005, 0.001)

numObs = 100
numVar = 5

#p = 1
#N = 1

for(p in 1:5){
    for(N in 1:100){
        u = rnorm(numObs, mean = 2, sd = 1)
        b = mvrnorm(numObs, mu = rep(0, numVar), Sigma = diag(0.1, numVar))
        x = u + b
        xsvd = svd(x)
        
        ## add noise
        temp = rep(NA, 5)
        temp[1] = xsvd$d[1]
        for(i in 2:numVar){
            temp[i] = xsvd$d[i] + rexp(1, rate = alpha[p]/xsvd$d[i])
        }
        xE = (xsvd$u %*% diag(temp) %*% t(xsvd$v))
        xEsvd = svd(xE)
        
        ## reverse noise
        temp2 = rep(NA, 5)
        temp2[1] = xEsvd$d[1]
        for(i in 2:numVar){
            temp2[i] = xEsvd$d[i] - rexp(1, rate = alpha[p]/xEsvd$d[i])
        }
        xRev = (xEsvd$u %*% diag(temp2) %*% t(xEsvd$v))
        
        xPC = PCA(x, graph = F)
        xRevPC = PCA(xRev, graph = F)
        xEPC = PCA(xE, graph = F)
        
        ####
        ## utility
        outUtil[N, p] = mean(abs(xPC$var$coord[,1] - xEPC$var$coord[,1]))
        outUtil[N + 100, p] = mean(abs(xPC$var$coord[,1] - xRevPC$var$coord[,1]))
        
        ## risk
        outRisk[N, p] = mean(rowMeans(abs(x - xE)))
        outRisk[N + 100, p] = mean(rowMeans(abs(x - xRev)))
        
        ##
        ##
        ## add noise
        temp = matrix(NA, nrow = numObs, ncol = numVar)
        temp[, 1] = xsvd$u[, 1]
        for(i in 2:numVar){
            temp[, i] = xsvd$u[, i] + runif(numObs, min = -epsilon[p], max = epsilon[p])
        }
        xE = (temp %*% diag(xsvd$d) %*% t(xsvd$v))
        
        xPC = PCA(x, graph = F)
        xEPC = PCA(xE, graph = F)
        
        ####
        ## utility
        outUtil[N + 200, p] = mean(abs(xPC$var$coord[,1] - xEPC$var$coord[,1]))
        
        ## risk
        outRisk[N + 200, p] = mean(rowMeans(abs(x - xE)))
        
        cat(p, N, '\n')
    }
}

riskHold = outRisk
utilHold = outUtil

combRisk = cbind(riskHold, outRisk[,1:5])
combUtil = cbind(utilHold, outUtil[,1:5])

#####
##### Summarize
#####
RUplot = data.frame(matrix(NA, ncol = 3, nrow = 45))
colnames(RUplot) = c("Risk", "Utility", "Type")
RUplot$Risk = c(scale(colMeans(combRisk[1:100,])), scale(colMeans(combRisk[101:200,])), 
                scale(colMeans(combRisk[201:300,])) )
RUplot$Utility = c(scale(colMeans(combUtil[1:100,])), scale(colMeans(combUtil[101:200,])), 
                   scale(colMeans(combUtil[201:300,])) )
RUplot$Type = rep(c("NoisyD", "ReverseD", "NoisyU"), each = 15)

RUplot$Risk = c(colMeans(combRisk[1:100,]), colMeans(combRisk[101:200,]), 
                colMeans(combRisk[201:300,]) )
RUplot$Utility = c(colMeans(combUtil[1:100,]), colMeans(combUtil[101:200,]), 
                   colMeans(combUtil[201:300,]) )

library(ggplot2)
p1 = ggplot(data = RUplot[c(1:13,31:43),], aes(x = -Risk, y = -Utility, color = Type)) + 
    geom_line() + geom_point()
p1

p2 = ggplot(data = RUplot[c(16:23,31:38),], aes(x = -Risk, y = -Utility, color = Type)) + 
    geom_line() + geom_point()
p2

p3 = ggplot(data = RUplot[c(1:8,16:23),], aes(x = -Risk, y = -Utility, color = Type)) + 
    geom_line() + geom_point()
p3
