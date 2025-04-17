# exercise 3
library(ggplot2)

#(a)
y <- c(0,2,5,5,7,1,4)
sigma <- 5

#the expression of posterior distribution unnormalized
posterior <- function(lambda,sigma,y){
  s <- sum(y)
  n <- length(y)
  res <- lambda^(s)*exp(-n*lambda-(lambda^2)/(2*sigma^2))
  return(res)
}

lambda_grid <- seq(0,7,by=0.05)

unnormalized <- posterior(lambda_grid,sigma,y)
normalized <- unnormalized/sum(unnormalized*0.05)

#plot(x=lambda_grid,y=normalized,type = "p")
df <- data.frame(lambda=lambda_grid,density=normalized)
ggplot(data=df,mapping = aes(x=lambda,y=density))+
  geom_point(color="#4E79A7")+
  labs(title = "Plot of posterior distribution ")+
  theme_classic()

#(b)
mode_index <- which.max(normalized)
mode <- lambda_grid[mode_index]
cat("The approximate posterior mode of lambda is",mode,"\n")



