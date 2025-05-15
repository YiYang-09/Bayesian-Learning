library("BayesLogit")
library("mvtnorm")
#(a)
#read data
data <- read.csv("Disease.csv")
data <- as.data.frame(data)

#normalize age, symtoms and white blood counts
age_mean <- mean(data$age)
age_sd <- sd(data$age)
duration_mean <- mean(data$duration_of_symptoms)
duration_sd <- sd(data$duration_of_symptoms)
white_mean <- mean(data$white_blood)
white_sd <- sd(data$white_blood)

data[,2] <- (data[,2]-age_mean)/age_sd
data[,4] <- (data[,4]-duration_mean)/duration_sd
data[,6] <- (data[,6]-white_mean)/white_sd

y <- data[,7] # response
X <- as.matrix(data[,1:6]) # Covariates


# prior sigma2
sigma2 <- 3^2

#number of iterations
num_iterations <- 1000


#Gibbs sampling
gibbs <- function(X,y,num_iterations,sigma2){
  if(nrow(X)!=length(y))
    stop("X and y should have the same length.")
  
  n <- nrow(X)
  p <- ncol(X)
  
  #initialize beta
  beta_samples <- matrix(NA,nrow = num_iterations,ncol = p)
  beta <- rep(0,p)
  
  
for (i in 1:num_iterations) {
  #latent variable omega
  #scale parameter for PG(1,x*beta)
  scalePar <- X%*%beta
  #construct latent variables for each observation thus num=n
  omega_variable <- rpg(num=n,h=rep(1,n),z=scalePar) #poly-Gamma sampling
  
  #update beta
  Omega_diagmatrix <- diag(omega_variable)
  
  V_inverse <- t(X) %*% Omega_diagmatrix %*%X +diag(1/sigma2,p)
  V <- solve(V_inverse) 
 
  k <- y-0.5
  b <- rep(0,p) #prior mean
  m <- V%*%(t(X)%*%k+solve(diag(1/sigma2,p))%*%b)           
  beta <- as.vector(rmvnorm(1,mean = m,sigma = V ))
  
  beta_samples[i,] <- beta
}
  return(beta_samples)
}

result <- gibbs(X=X,y=y,num_iterations=1000,sigma2 = 3^2)

#compute IFs ,lag.max = max_lag ,max_lag=100

Ifs <- function(chain){
  acf <- acf(chain,plot=F)$acf[-1]
  return(1+2*sum(acf))
}
  
IFs_a <- apply(result, 2,Ifs)

#compute Ifs using coda package
#library("coda")
#mcmc_result <- as.mcmc(result)
#ess <- effectiveSize(mcmc_result)
#ifactor <- nrow(result)/ess

#plot the trajectories of the sampled Markov chains

par(mfrow=c(2,3))
for (i in 1:ncol(result)) {
  plot(result[,i],type="l",col="#4E79A7",
       main=paste("Traceplot of beta",i),
       xlab=paste("beta",i),
       ylab="Iteartion",
       lwd=1)
}


#(b)
#m=10
result_10 <- gibbs(X=X[1:10,],y=y[1:10],num_iterations=1000,sigma2 = 3^2)
IFs_10 <- apply(result_10, 2, Ifs)

par(mfrow=c(2,3))
for (i in 1:ncol(result_10)) {
  plot(result_10[,i],type="l",col="#4E79A7",
       main=paste("Traceplot of beta",i),
       xlab=paste("beta",i),
       ylab="Iteartion",
       lwd=1)
}

#m=40

result_40 <- gibbs(X=X[1:40,],y=y[1:40],num_iterations=1000,sigma2 = 3^2)
IFs_40 <- apply(result_40, 2, Ifs)

par(mfrow=c(2,3))
for (i in 1:ncol(result_40)) {
  plot(result_40[,i],type="l",col="#4E79A7",
       main=paste("Traceplot of beta",i),
       xlab=paste("beta",i),
       ylab="Iteartion",
       lwd=1)
}

#m=80
result_80 <- gibbs(X=X[1:80,],y=y[1:80],num_iterations=1000,sigma2 = 3^2)
IFs_80 <- apply(result_80, 2, Ifs)

par(mfrow=c(2,3))
for (i in 1:ncol(result_80)) {
  plot(result_80[,i],type="l",col="#4E79A7",
       main=paste("Traceplot of beta",i),
       xlab=paste("beta",i),
       ylab="Iteartion",
       lwd=1)
}








