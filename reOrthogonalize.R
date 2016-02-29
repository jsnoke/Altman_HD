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