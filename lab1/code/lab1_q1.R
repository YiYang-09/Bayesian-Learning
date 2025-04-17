library(ggplot2)
#(a)
n <- 78
f <- 35
alpha0 <- 7
beta0 <- 7
alpha <- alpha0+(n-f)
beta <- beta0+f


set.seed(12345)
draws <- rbeta(n=10000,shape1 = alpha,shape2 = beta)

#derive the true mean and sd from beta distribution
true_mean <- alpha/(alpha+beta)
true_sd <- sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))

#derive the accumulating mean and sd values from random draws
accumulating_mean <- cumsum(draws)/(1:10000)
cum_sd <- function(data){
  n <- length(data)
  res <- numeric(n)
  for (i in 1:n) {
    res[i] <- sd(data[1:i])
  }
  return(res)
}
accumulating_sd <- cum_sd(draws)

#plot
#plot(accumulating_mean,type = "l",xlab = "number of draws",ylab = "mean",main = "Posterior Mean Convergence")
#abline(h=true_mean,col="red")

df_mean <- data.frame(x=1:length(accumulating_mean),y=accumulating_mean)
ggplot(data=df_mean,mapping = aes(x=x,y=y))+
  geom_line(color="#4E79A7")+
  geom_hline(yintercept = true_mean,color="red")
  labs(x="number of draws",y="mean",title = "Posterior Mean Convergence")+
  theme_minimal()
  
#(b)
set.seed(12345)
samples <- rbeta(10000,shape1 = alpha, shape2 = beta)  

# compute the posterior probability
prob <- mean(samples>0.5)

#derive the exact value from pbeta function
exactValue <- 1-pbeta(0.5,shape1 = alpha, shape2 = beta )
  
cat("The  posterior probability of theta > 0.5 :",prob,"\n")
cat("The exact value of theta > 0.5 from the Beta posterior :",exactValue,"\n")
  
#(c)
odds <- samples/(1-samples)
# set probability = TRUE(density histogram)
hist(odds,breaks = 50,probability = TRUE,main = "Posterior Distribution of Odds",
     xlab = "Odds",col ="#4E79A7" )
lines(density(odds),col="red",lwd=2)








