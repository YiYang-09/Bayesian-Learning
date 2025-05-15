#(a)

#AR(1)-process
set.seed(12345)
AR1 <- function(mu,phi,sigma2,t){
  points <- numeric(t)
  points[1] <- mu
  for (i in 2:t){
    xt <- mu+phi*(points[i-1]-mu)+rnorm(1,mean = 0,sd=sqrt(sigma2))
    points[i] <- xt
  }
return(points)
}

result1 <- AR1(mu=5,phi=0.8,sigma2 = 9,t=300)

phi_vals <- c(-0.9,-0.2,0,0.2,0.9)
par(mfrow=c(2,3))
set.seed(12345)
for (phi in phi_vals) {
  x <- AR1(mu=5,phi=phi,sigma2 = 9,t=300)
  plot(x,type="l",main=paste("phi = ",phi),
       xlab=" Time ",
       ylab=expression(x[t]),
       col="#4E79A7",
       lwd=1)
  
}

#(b)
library("rstan")

#phi=0.4
data1 <- AR1(mu=5,phi=0.4,sigma2=9,t=300)

stan_data1 <- list(
  T=length(data1),
  x=data1
)

stan_model_code <- "
data{
  int<lower=1> T; #time points
  vector[T] x;    #sequence
}

parameters {
  real mu;
  real<lower=-1, upper=1> phi;
  real<lower=0> sigma;  #! using sigma instead of sigma2, more stable
}

  
model{
  mu~normal(0,100);
  phi ~ uniform(-1, 1);
  sigma~cauchy(0,5);  #heavy tail
  x[1]~normal(mu,sigma);
  for (t in 2:T) {
    x[t]~ normal(mu+phi*(x[t-1]-mu),sigma);
  }
}  
"

fit1 <- stan(
  model_code = stan_model_code,
  data=stan_data1,
  iter=2000,
  chains=4,
  seed=12345
)



#phi=0.98
data2 <- AR1(mu=5,phi=0.98,sigma2=9,t=300)

stan_data2 <- list(
  T=length(data1),
  x=data2
)

fit2 <- stan(
  model_code = stan_model_code,
  data=stan_data2,
  iter=2000,
  chains=4,
  seed=12345
)

# mean of posterior, 95%ci,num of effective samples
sum_fit1 <-summary(fit1)$summary
fit1_df <- sum_fit1[1:3,c("mean","2.5%","97.5%","n_eff")]
print(fit1_df)

sum_fit2 <-summary(fit2)$summary
fit2_df <- sum_fit2[1:3,c("mean","2.5%","97.5%","n_eff")]
print(fit2_df)

#evaluate the convergence of tha samplers
print(fit1, pars=c("mu", "phi", "sigma"))
print(fit2, pars=c("mu", "phi", "sigma"))



#plot the joint posteriors of mu and phi
par(mfrow=c(1,1))
pairs(fit1, pars=c("mu", "phi"))
pairs(fit2, pars=c("mu", "phi"))

library(bayesplot)


posterior1 <- as.array(fit1)
mcmc_pairs(posterior1, pars = c("mu", "phi"))

  
posterior2 <- as.array(fit2)
mcmc_pairs(posterior2, pars = c("mu", "phi"))
 
