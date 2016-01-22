#----------------------------------------------------------------------------------------------------------------
prepared_data = read.delim("~/Box Sync/SDC Longitudinal Project/Towards Getting a Product/prepared_data.txt")

library(FactoMineR)
library(synthpop)

data(USArrests)
coo = PCA(USArrests,scale.unit=TRUE)


test = PCA(prepared_data, scale.unit = F, ncp = 10)

newPC = syn(test$ind$coord, method = "parametric", m = 5)

testSyn = test
testSyn$ind$coord = newPC$syn[[1]]

newVar = reconst(testSyn)
summary(newVar[, 1:5])
summary(prepared_data[, 1:5])

summary(reconst(test))

foo = prepared_data[, c(1:3, 5, 11, 13, 15:16, 18, 20)]
summary(foo)
test = PCA(foo, scale.unit = T)

newFoo = reconst(test, ncp = 3)
summary(newFoo)


#----------------------------------------------------------------------------------------------------------------
uDat = rnorm(100, mean = 5, sd = 3)
zDat = mvrnorm(n = 100, mu = rep(0, 3), Sigma = diag(2, 3))
xDat = uDat + zDat
colnames(xDat) = paste("x", 1:3, sep = "")
summary(xDat)

testPCA = PCA(xDat, scale.unit = T)
otherTest = prcomp(xDat)

newData = reconst(testPCA, ncp = 1)
summary(newData)
summary(testData)


##mult Imp
synthX = vector("list", 5)

for(j in 1:5){
    synthU = mean(xDat)
    synthZ = matrix(NA, nrow = 100, ncol = 3)
    for(i in 1:3){
        synthZ[, i] = rnorm(100, mean = mean(xDat[, i]) - synthU, sd = sd(xDat[, i]))
    }
    synthX[[j]] = synthU + synthZ
    colnames(synthX[[j]]) = paste("x", 1:3, sep = "")
}

summary(synthX)

foo

#--------------------------------------------------
x = mvrnorm(n = 1000, mu = c(25, 0, 0, 0), Sigma = matrix(c(10, 0, 0, 0,
                                                            0, 0, 0, 0,
                                                            0, 0, 0, 0,
                                                            0, 0, 0, 0), nrow = 4, ncol = 4))
summary(x)

#beta = matrix(c(2, 0, 0, 0))
epsilon = mvrnorm(n = 1000, mu = c(0, 0, 0, 0), Sigma = matrix(c(1, 0, 0, 0,
                                                                 0, 1, 0, 0,
                                                                 0, 0, 1, 0,
                                                                 0, 0, 0, 1), nrow = 4, ncol = 4))

w = x + epsilon
summary(w) ## true mean

standCI = matrix(NA, nrow = 100, ncol = 3)
PCCI = matrix(NA, nrow = 100, ncol = 3)
for(i in 1:100){
    sampW = w[sample(1:1000, 100),]
    summary(sampW)
    
    ##Standard Way
    newWish = rWishart(1, 99, solve(cov(sampW))/99 )
    newMu = mvrnorm(n = 1, apply(sampW, 2, mean), solve(newWish[,,1])/100)
    newW = mvrnorm(n = 1000, newMu, solve(newWish[,,1]))
    
    sampNewW = newW[sample(1:1000, 100),]
    #summary(sampNewW)
    standCI[i,] = c(mean(sampNewW[,1]), mean(sampNewW[,1]) - 1.96*sd(sampNewW[,1]), mean(sampNewW[,1]) + 1.96*sd(sampNewW[,1]))
    
    
    ##PCA Way
    wPC = PCA(sampW, scale.unit = F)
    #summary(wPC$ind$coord)
    
    pcWish = rWishart(1, 99, solve(cov(wPC$ind$coord)[1,1])/99 )
    pcMu = mvrnorm(n = 1, mean(wPC$ind$coord[,1]), solve(pcWish[,,1])/100)
    newPC = mvrnorm(n = 1000, pcMu, solve(pcWish[,,1]))
    
    sampNewPC = newPC[sample(1:1000, 100),]
    
    newWPC = wPC
    newWPC$ind$coord[, 1] = sampNewPC
    #summary(reconst(newWPC))
    
    foobar = reconst(newWPC, ncp = 1)[,1]
    PCCI[i,] = c(mean(foobar), mean(foobar) - 1.96*sd(foobar), mean(foobar) + 1.96*sd(foobar))
    
    
}


#-----------------
W1_public_20130621 = read.csv("~/Box Sync/SDC Longitudinal Project/CSV DATA FILES/W1_public_20130621.csv", comment.char="#")

W2_public_111103 = read.csv("~/Box Sync/SDC Longitudinal Project/CSV DATA FILES/W2_public_111103.csv", comment.char="#")


View(W1_public_20130621)

foo = sapply(W1_public_20130621, class)

foo

numNumeric = sapply(W1_public_20130621, is.numeric) ## 523 var., mostly ordinal and categorical
numFactor = sapply(W1_public_20130621, is.factor) ## 742 var., categorical, few ordinal
numLogical = sapply(W1_public_20130621, is.logical) ## 127 var., all NAs

dim(W2_public_111103)

foo = sapply(W2_public_111103, class)

foo

numNumeric = sapply(W2_public_111103, is.numeric) ## 419 var., mostly ordinal and categorical
numFactor = sapply(W2_public_111103, is.factor) ## 317 var., categorical, few ordinal
numLogical = sapply(W2_public_111103, is.logical) ## 93 var., all NAs


















