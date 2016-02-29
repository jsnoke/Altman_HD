#######################

##Risk measure - DP
##Following Jiang et. al 2015

#######################
library(FactoMineR)
library(VGAM)

numObs = 100
numVar = 1000
beta = c(0.5, 3)

x = rnorm(numObs, mean = 2, sd = 1)
y = matrix(NA, nrow = numObs, ncol = numVar)
for(i in 1:ncol(y)){
    y[, i] = beta[1] + beta[2]*x + rnorm(numObs)
}

d = c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90,
      99, 100, 101, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
resultsR = vector("list", length(d))
resultsE = vector("list", length(d))
for(a in 1:length(d)){
    DF = y[ , 1:d[a]]
    DFE = DF
    
    epsilon = seq(0.01, 1, 0.09)
    resultsR[[a]] = matrix(NA, nrow = d[a], ncol = length(epsilon))
    resultsE[[a]] = matrix(NA, nrow = d[a], ncol = length(epsilon))
    for(b in 1:length(epsilon)){
        for(i in 1:ncol(DF)){
            DFE[, i] = DF[, i] + rlaplace(numObs, 0, (3*d[a]) / epsilon[b] )
        }
        boo = PCA(DF, graph = F)
        foo = PCA(DFE, graph = F)
        
        resultsR[[a]][, b] = dimdesc(boo, axes = 1, proba = 0.05)$Dim.1$quanti[, 1]
        resultsE[[a]][, b] = dimdesc(foo, axes = 1, proba = 1)$Dim.1$quanti[, 1]
        
        cat(b)
    }
    cat(a, "\n")
}

par(mfcol = c(2, 5))
for(c in 1:5){
    boxplot(c(resultsR[[c]]))
    boxplot(c(resultsE[[c]]))
}




#factanal(DF, 1)
#factanal(DFE, 1)


#####
##
#####
epsilon = 0.1
noiseL = matrix(NA, nrow = d, ncol = d)
noiseL[upper.tri(noiseL, diag = T)] = rlaplace( ((d^2 + d) / 2), 0, (2 * d) / (numObs*epsilon) )
for(i in 1:d){
    for(j in 1:d){
        noiseL[j, i] = noiseL[i, j]
    }
}

noiseCOV = (t(XDF) %*% XDF) / numObs + noiseL
noiseDecomp = eigen(noiseCOV)
noiseX = origSVD$u %*% diag(noiseDecomp$values) %*% t(noiseDecomp$vectors)

#eigen(cov(XDF))
#eigen(noiseCOV)

origPCA = princomp(XDF)
cor(origPCA$scores[,1], x1)



## wishart method
noiseCOV2 = rWishart(1, d, C)[, , 1] + cov(XDF)
noisePCA = princomp(covmat = noiseCOV2)


##### make orthogonal
setEign = 3 / (2 * numObs * epsilon)

tmp <- rnorm(2) 
tmp.qr <- qr(tmp) 
tmp.complete <- qr.Q(tmp.qr, complete=TRUE) 

C = tmp.complete %*% diag(x = setEign, 2) %*% t(tmp.complete)


#############
scaleData = scale(XDF)
origSVD = svd(scaleData)
gsU = c(sqrt(numObs) / abs((origSVD$d[1] - origSVD$d[2])))
gsLam = sqrt(2*ncol(scaleData))
epsilon = 1
noiseU = origSVD$u[, 2:3, drop = F] + rlaplace(numObs, scale = gsU / epsilon )
noiseLam = origSVD$d[2:3] + rlaplace(1, scale = gsLam / epsilon )

tmp = qr(noiseU)
orthU = qr.Q(tmp, complete=TRUE)[,1:2]

summary(orthU %*% diag(noiseLam, 2) %*% t(origSVD$v[, 2:3, drop = F]))
summary(origSVD$u[, 2:3] %*% diag(origSVD$d[2:3], 2) %*% t(origSVD$v[, 2:3]))
plot(orthU %*% diag(noiseLam, 2) %*% t(origSVD$v[, 2:3, drop = F]))
plot(origSVD$u[, 2:3] %*% diag(origSVD$d[2:3], 2) %*% t(origSVD$v[, 2:3]))
cor(orthU %*% diag(noiseLam, 2) %*% t(origSVD$v[, 2:3, drop = F]))
cor(origSVD$u[, 2:3] %*% diag(origSVD$d[2:3], 2) %*% t(origSVD$v[, 2:3]))




















