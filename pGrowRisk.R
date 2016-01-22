library(MASS)
library(synthpop)
library(FactoMineR)
#########
## new risk eval
#########

rowDiff = function(df, row){
    return(abs(row - df))
}

closest = function(list, tol){
    return(which(rowSums((list < tol)) == max(rowSums((list < tol)))))
}

numSim = 100
numObs = 50
numVar = c(5, 25, 51, 500)
matches = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
               "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matchMSE = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                 "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matches2 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
               "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matchMSE2 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matches3 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matchMSE3 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                 "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matches4 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matchMSE4 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                 "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matches5 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))
matchMSE5 = list("five" = vector("list", numSim), "twentyfive" = vector("list", numSim),
                 "fiftyone" = vector("list", numSim), "hundred" = vector("list", numSim))


for(i in 1:length(numVar)){
    for(j in 1:numSim){
        u = rnorm(numObs, mean = 2, sd = 1)
        z = mvrnorm(numObs, mu = rep(0, numVar[i]), Sigma = diag(0.1, numVar[i]))
        x = u + z
        summary(x)
        colnames(x) = paste("X", 1:numVar[i], sep = "")
        
        temp = vector("list", 5)
        for(b in 1:5){
            synPC = PCA(x, scale.unit = F, ncp = 1)
            temp[[b]] = reconst(synPC, ncp = 1) 
            #for(c in 1:numVar[i]){
            #    temp[[b]][, c] = temp[[b]][, c] + rnorm(numObs, 0, 2)
            #}
        }
        #temp = syn(x, m = 5)$syn
        synX = rbind(temp[[1]], temp[[2]], temp[[3]], temp[[4]], temp[[5]])
        
        #foo = t(apply(as.matrix(synX$syn[, 1:(numVar[i] - 1)]), 1, rowDiff, row = x[1, 1:(numVar[i] - 1), drop = F]))
        foo = t(apply(as.matrix(synX[, 1:(numVar[i] - 1)]), 1, rowDiff, row = x[1, 1:(numVar[i] - 1), drop = F]))
        
        #matches3[[i]][[j]] = rowSums(foo < sd(x))
        
        closeInd = closest(foo, sd(x))
        matches5[[i]][[j]] = as.numeric(rowSums(foo < sd(x))[closeInd])[1]
        
        #allMatch =  abs(synX$syn[closeInd, numVar[i]] - x[1, numVar[i]])
        allMatch =  abs(synX[closeInd, numVar[i]] - x[1, numVar[i]])
        
        if(is.na(var(allMatch)) == TRUE){
            matchMSE5[[i]][[j]] = mean(allMatch)^2
        } else
            matchMSE5[[i]][[j]] = mean(allMatch)^2 + var(allMatch)
        #hist(matchDist)        
    }
}

#####
summary(unlist(matches$five))
summary(unlist(matches$twentyfive))
summary(unlist(matches$fiftyone))
summary(unlist(matches$hundred))

summary(unlist(matchMSE$five))
summary(unlist(matchMSE$twentyfive))
summary(unlist(matchMSE$fiftyone))
summary(unlist(matchMSE$hundred))

##
summary(unlist(matches2$five))
summary(unlist(matches2$twentyfive))
summary(unlist(matches2$fiftyone))
summary(unlist(matches2$hundred))

summary(unlist(matchMSE2$five))
summary(unlist(matchMSE2$twentyfive))
summary(unlist(matchMSE2$fiftyone))
summary(unlist(matchMSE2$hundred))

##
summary(unlist(matches3$five))
summary(unlist(matches3$twentyfive))
summary(unlist(matches3$fiftyone))
summary(unlist(matches3$hundred))

summary(unlist(matchMSE3$five))
summary(unlist(matchMSE3$twentyfive))
summary(unlist(matchMSE3$fiftyone))
summary(unlist(matchMSE3$hundred))

##
summary(unlist(matches5$five))
summary(unlist(matches5$twentyfive))
summary(unlist(matches5$fiftyone))
summary(unlist(matches5$hundred))

summary(unlist(matchMSE5$five))
summary(unlist(matchMSE5$twentyfive))
summary(unlist(matchMSE5$fiftyone))
summary(unlist(matchMSE5$hundred))

##
par(mfrow = c(2,2))
plot(unlist(matches3$five), unlist(matchMSE3$five))
plot(unlist(matches3$twentyfive), unlist(matchMSE3$twentyfive))
plot(unlist(matches3$fiftyone), unlist(matchMSE3$fiftyone))
plot(unlist(matches3$hundred), unlist(matchMSE3$hundred))

cor(unlist(matches3$five), unlist(matchMSE3$five))
cor(unlist(matches3$twentyfive), unlist(matchMSE3$twentyfive))
cor(unlist(matches3$fiftyone), unlist(matchMSE3$fiftyone))
cor(unlist(matches3$hundred), unlist(matchMSE3$hundred))

##
par(mfrow = c(1,3))
plot(c(unlist(matches$five), unlist(matches$twentyfive), unlist(matches$fiftyone), unlist(matches$hundred)),
    c(unlist(matchMSE$five), unlist(matchMSE$twentyfive), unlist(matchMSE$fiftyone), unlist(matchMSE$hundred)))
plot(c(unlist(matches3$five), unlist(matches3$twentyfive), unlist(matches3$fiftyone), unlist(matches3$hundred)),
       c(unlist(matchMSE3$five), unlist(matchMSE3$twentyfive), unlist(matchMSE3$fiftyone), unlist(matchMSE3$hundred)))
plot(c(unlist(matches4$five), unlist(matches4$twentyfive), unlist(matches4$fiftyone), unlist(matches4$hundred)),
     c(unlist(matchMSE4$five), unlist(matchMSE4$twentyfive), unlist(matchMSE4$fiftyone), unlist(matchMSE4$hundred)))


cor(c(unlist(matches3$five), unlist(matches3$twentyfive), unlist(matches3$fiftyone), unlist(matches3$hundred)),
     c(unlist(matchMSE3$five), unlist(matchMSE3$twentyfive), unlist(matchMSE3$fiftyone), unlist(matchMSE3$hundred)))
cor(c(unlist(matches4$five), unlist(matches4$twentyfive), unlist(matches4$fiftyone), unlist(matches4$hundred)),
     c(unlist(matchMSE4$five), unlist(matchMSE4$twentyfive), unlist(matchMSE4$fiftyone), unlist(matchMSE4$hundred)))

cor(c(unlist(matches$five), unlist(matches$twentyfive), unlist(matches$fiftyone), unlist(matches$hundred)),
     c(unlist(matchMSE$five), unlist(matchMSE$twentyfive), unlist(matchMSE$fiftyone), unlist(matchMSE$hundred)))



#####
#####
#####
numVar = 5
#numObs = 100
u = rnorm(numObs, mean = 2, sd = 1)
z = mvrnorm(numObs, mu = rep(0, numVar), Sigma = diag(0.1, numVar))
x = u + z

foo = svd(x)

summary(foo$d[-1])
summary(rowMeans(foo$u[, -1]))
summary(rowMeans(foo$v[, -1]))
hist(foo$d[-1])

foo$d
sqrt(foo$d)
sum(foo$d)



(sqrt(foo$d) - mean(sqrt(foo$d))) / sd(sqrt(foo$d))
head(foo$u)
foo$v

blah = princomp(x)
blah$sdev
blah$loadings
plot(blah)

head(blah$scores)
head(foo$u %*% diag(foo$d))

eH = foo$u[, 2:(numVar)] %*% diag(foo$d[2:numVar]) %*% t(foo$v[, 2:numVar])
eI = foo$u[, 1] %*% diag(foo$d[1], nrow = 1) %*% t(foo$v[, 1])
#head((eI + eH)[, 1:5])
head(x[, 1:5])
head(eI[, 1:5])
head(eH[, 1:5])

summary(rowMeans(eH)^2)

cov(eI)
cov(eH)

#####
##### Altman suggest - 11.19.15
#####
numVar = 500
numObs = 100
u = rnorm(numObs, mean = 2, sd = 1)
b = mvrnorm(numObs, mu = rep(0, numVar), Sigma = diag(0.1, numVar))
x = u + b

#sig = diag(1, numVar) + rnorm(numVar^2, 0, 0.1)
#x = mvrnorm(numObs, mu = rep(0, numVar), Sigma = sig) + rnorm(numObs, 0, 5)
summary(x)
cov(x)
cor(x)

xPca = PCA(x, ncp = numVar, scale.unit = F)
xPca$eig

head(xPca$svd$U%*%diag(xPca$svd$vs)%*%t(xPca$svd$V))
head(x)

head(xPca$svd$V %*% diag(xPca$svd$vs^2) %*% t(xPca$svd$V))
head(t(x) %*% x) / numObs
#head(t(xPca$var$coord) %*% diag(xPca$eig$eigenvalue) %*% xPca$var$coord)

z = matrix(NA, nrow = numObs, ncol = numVar)
for(i in 1:numObs){
    for(j in 1:numVar){
        z[i, j] = rnorm(1, 0, 1)
    }
}
summary(z)

summary(z %*% diag(xPca$svd$vs) %*% t(xPca$svd$V))
summary(x)


