#####
## risk plots
#####
load("~/Documents/Altman_HD/R/miNoise4.RData")
load("~/Documents/Altman_HD/R/partNoise4.RData")
load("~/Documents/Altman_HD/R/noErr4.RData")

library(ggplot2)
library(grid)

plotDF = data.frame(matrix(NA, nrow = 3*10*6, ncol = 5))
colnames(plotDF) = c("Factor1", "Factor2", "Noise", "p", "epsilon")

plotDF[, "Factor1"] = c(resultsMu1, results2Mu1, results3Mu1)
plotDF[, "Factor2"] = c(resultsMu2, results2Mu2, results3Mu2)
plotDF[, "Noise"] = rep(c("Single Impute", "Noise Only", "Multiple Imputes"), each = 60)
plotDF[, "p"] = rep(rep(p, each = 10), 3)
plotDF[, "epsilon"] = rep(rep(epsilon, 6), 3)


p1 = ggplot(data = plotDF[plotDF$p == 25, ], aes(x = epsilon, y = Factor1, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor1)) + ggtitle("p = 25") + 
    ylab("Mean Absolute Loadings Difference") + theme_bw()
p2 = ggplot(data = plotDF[plotDF$p == 101, ], aes(x = epsilon, y = Factor1, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor1)) + ggtitle("p = 101") + 
    ylab("Mean Absolute Loadings Difference") + theme_bw()
p3 = ggplot(data = plotDF[plotDF$p == 500, ], aes(x = epsilon, y = Factor1, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor1)) + ggtitle("p = 500") + 
    ylab("Mean Absolute Loadings Difference") + theme_bw()
p4 = ggplot(data = plotDF[plotDF$p == 1000, ], aes(x = epsilon, y = Factor1, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor1)) + ggtitle("p = 1000") + 
    ylab("Mean Absolute Loadings Difference") + theme_bw()

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1) )
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2) )
print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1) )
print(p4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2) )


pp1 = ggplot(data = plotDF[plotDF$p == 25, ], aes(x = epsilon, y = Factor2, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor2)) + ggtitle("p = 25") + 
    ylab("Mean Absolute Loadings Difference")
pp2 = ggplot(data = plotDF[plotDF$p == 101, ], aes(x = epsilon, y = Factor2, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor2)) + ggtitle("p = 101") + 
    ylab("Mean Absolute Loadings Difference")
pp3 = ggplot(data = plotDF[plotDF$p == 500, ], aes(x = epsilon, y = Factor2, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor2)) + ggtitle("p = 500") + 
    ylab("Mean Absolute Loadings Difference")
pp4 = ggplot(data = plotDF[plotDF$p == 1000, ], aes(x = epsilon, y = Factor2, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor2)) + ggtitle("p = 1000") + 
    ylab("Mean Absolute Loadings Difference")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(pp1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1) )
print(pp2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2) )
print(pp3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1) )
print(pp4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2) )


#####
## syn plots
#####

load("~/Documents/Altman_HD/R/noHigh.RData")

p = c(5, 10, 25, 50, 75, 90)
resultsMu1 = matrix(NA, ncol = 6, nrow = 2*1000)
resultsMu2 = matrix(NA, ncol = 6, nrow = 2*1000)
#resultsSd1 = matrix(NA, ncol = 6, nrow = 2*1000)
#resultsSd2 = matrix(NA, ncol = 6, nrow = 2*1000)
for(c in 1:2){
    for(b in 1:6){
        #resultsMu1[c, b] = 0
        #resultsMu2[c, b] = 0
        #resultsSd1[c, b] = 0
        #resultsSd2[c, b] = 0
        for(a in 1:1000){
            hold = test[a, b][[1]]
            resultsMu1[(a + (c - 1) * 1000), b] = mean(abs(hold[, 1] - hold[, (1 + c * 2)]))
            #resultsMu1[c, b] = (mean(abs(hold[, 1] - hold[, (1 + c * 2)])) + resultsMu1[c, b] * (a - 1)) / a
            resultsMu2[(a + (c - 1) * 1000), b] = mean(abs(hold[, 2] - hold[, (2 + c * 2)]))
            #resultsMu2[c, b] = (mean(abs(hold[, 2] - hold[, (2 + c * 2)])) + resultsMu2[c, b] * (a - 1)) / a
            
            #resultsSd1[c, b] = (sd(abs(hold[, 1] - hold[, (1 + c * 2)])) + resultsSd1[c, b] * (a - 1)) / a
            
            #resultsSd2[c, b] = (sd(abs(hold[, 2] - hold[, (2 + c * 2)])) + resultsSd2[c, b] * (a - 1)) / a
        }
    }
}

library(ggplot2)
library(grid)

plotDF = data.frame(matrix(NA, nrow = 2000*6, ncol = 4))
colnames(plotDF) = c("Factor1", "Factor2", "Synthesis", "p")

plotDF[, "Factor1"] = c(resultsMu1)
plotDF[, "Factor2"] = c(resultsMu2)
plotDF[, "Synthesis"] = rep(rep(c("Sequential Synthesis", "SVD Synthesis"), each = 1000), 6)
plotDF[, "p"] = rep(p, each = 2000)


p1 = ggplot(data = plotDF, aes(x = p, y = Factor1, color = Synthesis)) + ggtitle("Utility - First Factor") + 
    ylab("Mean Absolute Loadings Difference") + stat_smooth(method = "lm", se = F) + theme_bw()
#p1
p2 = ggplot(data = plotDF, aes(x = p, y = Factor2, color = Synthesis)) +  ggtitle("Utility - Second Factor") + 
    ylab("Mean Absolute Loadings Difference") + stat_smooth(method = "lm", se = F) + theme_bw()
#p2

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1) )
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2) )

#####
## lowD risk plots
#####

load("~/Documents/Altman_HD/R/nonHighRisk.RData")

library(ggplot2)
library(grid)

p = c(5, 25, 50, 101, 500, 1000)

plotDF = data.frame(matrix(NA, nrow = 6*3, ncol = 4))
colnames(plotDF) = c("Total", "Prob", "p", "Neighborhood.Size")

plotDF[, "Total"] = c(colMeans(nonHighRisk[[1]]), colMeans(nonHighRisk[[2]]), colMeans(nonHighRisk[[3]]))
plotDF[, "Prob"] = c(colMeans(nonHighRisk[[4]]), colMeans(nonHighRisk[[5]]), colMeans(nonHighRisk[[6]]))
plotDF[, "p"] = rep(p, 3)
plotDF[, "Neighborhood.Size"] = rep(c("0.01", "0.1", "0.05"), each = 6)


p3 = ggplot(data = plotDF, aes(x = p, y = Total, color = Neighborhood.Size)) + 
    ggtitle("Risk - Individual Matches") + ylim(range(plotDF$Total)) + geom_line() +
    #stat_smooth(method = "lm", se = F) +
    ylab("Max Var. Matches for a Released Row") + theme_bw()
#p1

p4 = ggplot(data = plotDF, aes(x = p, y = Prob, color = Neighborhood.Size)) + 
    ggtitle("Risk - Individual Matches") + ylim(range(plotDF$Prob)) + geom_line() +
    #stat_smooth(method = "lm", se = F) +
    ylab("Probability of Unique Maximum Match") + theme_bw()

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1) )
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2) )
print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1) )
print(p4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2) )

#####
## theoretical ru map
#####

plotDF = data.frame(matrix(NA, nrow = 30, ncol = 3))
colnames(plotDF) = c("Risk", "Utility", "Theoretical.SDC.Method")

plotDF[, "Risk"] = c(c(0, 0, 0, 0.1, 0.1, 0.2, 0.2, 0.4, 0.6, 0.9),
                     c(0, 0.2, 0.4, 0.6, 0.8, 0.8, 0.9, 1, 1, 1),
                     c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
plotDF[, "Utility"] = c(c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 0.85, 0.9, 0.95),
                     c(0, 0.2, 0.4, 0.6, 0.8, 0.8, 0.9, 1, 1, 1),
                     c(0, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.2))
plotDF[, "Theoretical.SDC.Method"] = factor(c(rep(1:3, each = 10)))

ggplot(data = plotDF, aes(y = Utility, x = Risk, color = Theoretical.SDC.Method)) + stat_smooth(se = F) + 
    theme_bw() + theme(axis.title = element_text(size=20)) + ylab("Utility") + xlab("Risk")







