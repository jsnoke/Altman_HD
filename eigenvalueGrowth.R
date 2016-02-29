#######################

##
##

#######################
library(FactoMineR)
library(VGAM)

summary(rlaplace(500, 0, 0.5))

numObs = 100
numVar = 1000
beta = c(0.5, 3)

x = rnorm(numObs, mean = 2, sd = 1)
y = matrix(NA, nrow = numObs, ncol = numVar)
for(i in 1:ncol(y)){
    y[, i] = beta[1] + beta[2]*x + rnorm(numObs)
}

DF = scale(cbind(x, y[ , 1:1000]))

foo = svd(DF)
sum(foo$d[-1]^2)
sum(foo$d[-1]^2) / sum(foo$d^2)

boo = PCA(DF)
