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
    ylab("Mean Absolute Loadings Difference")
p2 = ggplot(data = plotDF[plotDF$p == 101, ], aes(x = epsilon, y = Factor1, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor1)) + ggtitle("p = 101") + 
    ylab("Mean Absolute Loadings Difference")
p3 = ggplot(data = plotDF[plotDF$p == 500, ], aes(x = epsilon, y = Factor1, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor1)) + ggtitle("p = 500") + 
    ylab("Mean Absolute Loadings Difference")
p4 = ggplot(data = plotDF[plotDF$p == 1000, ], aes(x = epsilon, y = Factor1, color = Noise)) + 
    geom_line() + ylim(range(plotDF$Factor1)) + ggtitle("p = 1000") + 
    ylab("Mean Absolute Loadings Difference")

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





