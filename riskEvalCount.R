#########
## new risk eval
#########
library(MASS)

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

synX = syn(xDat1, m = 1, method = "parametric")

##priors
meanY = mean(synX$syn[, 1])
foo = glm(X1 ~ ., data = synX$syn[-1,], family = poisson)
dpois(round(predict(foo, type = "response", newdata = synX$syn[1,])), lambda = meanY)












