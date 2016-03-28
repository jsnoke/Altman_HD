library(MASS)
x = mvrnorm(100, mu = rep(0, 4), Sigma = matrix(c(1, rep(-0.2, 3),
                                                  rep(-0.2, 1), 1, rep(-0.2, 2),
                                                  rep(-0.2, 2), 1, rep(-0.2, 1),
                                                  rep(-0.2, 3), 1), nrow = 4, ncol = 4, byrow = T) )
x = scale(x)
t(x)%*%x
c = solve(t(x)%*%x)

u = x%*%sqrt(c)
summary(u)
t(u)%*%u

tmp = qr(noiseU)
orthU = qr.Q(tmp, complete=FALSE)

#####
###
#####
library(FactoMineR)
library(VGAM)
library(MASS)

numObs = 100
numVar = 1000
beta = c(0.5, 3)

x = mvrnorm(numObs, mu = rep(2, 2), Sigma = matrix(c(1, 0, 0, 1), nrow = 2) )
y = matrix(NA, nrow = numObs, ncol = numVar)
for(i in 1:ncol(y)){
    if(runif(1) > 0.5)
        y[, i] = beta[1] + beta[2]*x[, 1] + rnorm(numObs)
    else
        y[, i] = beta[2] + beta[1]*x[, 2] + rnorm(numObs)
}

d = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90,
      99, 100, 101, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
DF = y[ , 1:d[a]]
#DFE = DF

epsilon = seq(0.01, 1, 0.09)

scaleData = scale(DF)
origSVD = svd(DF)
origSVD$d
#a = 1
gsU = c(sqrt(numObs) * 3 * sqrt(d[a]) / (origSVD$d[1] - origSVD$d[2]))
gsLam = sqrt(1) * 3 * sqrt(d[a])

noiseU = matrix(NA, nrow = numObs, ncol = 2)
noiseLam = rep(NA, 2)
for(i in 1:2){
    noiseU[, i] = scale(origSVD$u[, i, drop = F] + rlaplace(numObs, scale = gsU / epsilon[b] ))
    noiseLam[i] = origSVD$d[i] + rlaplace(1, scale = gsLam / epsilon[b] )
}

DFE = orthU %*% diag(noiseLam, 2) %*% origSVD$v[1:2, ]
test = origSVD$u[, 1:2] %*% diag(origSVD$d[1:2], 2) %*% origSVD$v[1:2, ]

boo = PCA(DF, graph = F)
foo = PCA(DFE, graph = F)

loadsR = dimdesc(boo, axes = 1, proba = 1)$Dim.1$quanti[, 1]
loadsE = dimdesc(foo, axes = 1, proba = 1)$Dim.1$quanti[, 1]




