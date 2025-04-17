library(bayestestR)
#(a)
y <- c(22,33,31,49,65,17,24)
mu <- 3.65
n <- length(y)
#calculate tau square
tauSquare <- sum((log(y)-mu)^2)/n

set.seed(12345)

samples <- (n*tauSquare)/rchisq(10000,df=n)
hist(samples,breaks = 50,main = "Posterior of sigma suqare",xlab = expression(sigma^2))

#(b)
sigma_samples <- sqrt(samples)
G <- 2*pnorm(sigma_samples/(sqrt(2)))-1

#(c)
mean_G <- mean(G)
sd_G <- sd(G)
lower_bound <- mean_G-1.96*sd_G
upper_bound <- mean_G+1.96*sd_G
cat("The Bayesian 95% cedible interval for G is [", lower_bound,upper_bound,"]","\n")

#(d)
HPDI <- hdi(G,ci=0.95)
cat("The Bayesian 95% highest posterior density interval for G is [",HPDI[[2]],HPDI[[3]],"]","\n")





