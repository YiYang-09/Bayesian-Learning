# Computer Lab2 - 732A73
# Helena Llorens Lluís (hllor282), Yi Yang (yiyan338)
#---------------------------------------------------------------------------------------------------------

# ASSIGNMENT 1

# load data
data <- read.csv("D:\\liu\\BL\\Bayesian-Learning\\lab2\\temp_linkoping.csv")
n <- nrow(data)
y <- data$temp
x <- data$time

# (a)
# prior hyperparameters
mu0 <- t(c(20, -100, 100))
omega0 <- 50 * diag(3)
nu0 <- 1
sigma0 <- 1

# scale inverse chi function
scale_inv_chi <- function(nu0, sigma0, N) {
  sigma0 <- as.numeric(sigma0)
  samples <- (nu0 * sigma0) / rchisq(N, df = nu0)
  return(samples)
}

# prior distributions
N <- 50
set.seed(123)
prior_sigma <- scale_inv_chi(nu0, sigma0, N)
prior_beta <- matrix(NA, nrow = N, ncol = 3)
for(i in 1:N){
  prior_beta[i, ] <- mvtnorm::rmvnorm(n = 1, mean = mu0, sigma = prior_sigma[i] * solve(omega0))
}

# Plot prior regression curves
plot(data$time, y, col = "#BAB0AC", pch = 20, 
     xlab = "Time", ylab = "Temperature", ylim = c(-20, 30))
for (i in 1:N) {
  curve <- prior_beta[i, 1] + prior_beta[i, 2] * data$time + prior_beta[i, 3] * data$time^2
  lines(data$time, curve, col = "#4E79A7")
}


# (b)
# posterior distribution function
X <- matrix(c(rep(1, n), data$time, data$time^2), ncol = 3)
Y <- y
posterior <- function(mu0, omega0, nu0, sigma0, X, Y, N, p){
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  mu_n <- solve(t(X) %*% X + omega0) %*% (t(X) %*% X %*% beta + omega0 %*% matrix(mu0, ncol = 1))
  omega_n <- t(X) %*% X + omega0
  nu_n <- nu0 + nrow(X)
  sigma_n <- (nu0 * sigma0 + 
                (t(Y) %*% Y + t(matrix(mu0, ncol = 1)) %*% omega0 %*% matrix(mu0, ncol = 1) - 
                   t(mu_n) %*% omega_n %*% mu_n))/nu_n
  post_sigma <- scale_inv_chi(nu_n, sigma_n, N)
  post_beta <- matrix(NA, nrow = N, ncol = p + 1)
  for(i in 1:N){
    post_beta[i, ] <- mvtnorm::rmvnorm(n = 1, mean = mu_n, sigma = post_sigma[i] * solve(omega_n))
  }
  
  return(list(sigma = post_sigma, beta = post_beta))
}

set.seed(123)
post <- posterior(mu0, omega0, nu0, sigma0, X, Y, N, 2)

# (i) Histograms of the magrinal posterior for each parameter
hist(post$sigma, col = "#4E79A7", xlab = expression(sigma), main = "")
hist(post$beta[, 1], col = "#4E79A7", xlab = expression(beta[0]), main = "")
hist(post$beta[, 2], col = "#4E79A7", xlab = expression(beta[1]), main = "")
hist(post$beta[, 3], col = "#4E79A7", xlab = expression(beta[2]), main = "")

# (ii)
# draws from the function
time_grid <- seq(min(data$time), max(data$time), length.out = 100)
f_pred_draws <- matrix(NA, nrow = N, ncol = length(time_grid))
for (i in 1:N) {
  f_mean <- post$beta[i, 1] + post$beta[i, 2] * time_grid + post$beta[i, 3] * time_grid^2
  f_pred_draws[i, ] <- rnorm(length(time_grid), mean = f_mean, sd = sqrt(post$sigma[i]))
}

# compute median and CI
pred_median <- apply(f_pred_draws, 2, median)
pred_lb <- apply(f_pred_draws, 2, quantile, probs = 0.05)
pred_ub <- apply(f_pred_draws, 2, quantile, probs = 0.95)

# Plot data
plot(x, y, col = "#BAB0AC", pch = 20,
     main = "",
     xlab = "Time", ylab = "Temperature")

# Add posterior median curve
lines(time_grid, pred_median, col = "#4E79A7", lwd = 2)

# Add 90% posterior credible interval
lines(time_grid, pred_lb, col = "#E15759", lty = 2)
lines(time_grid, pred_ub, col = "#E15759", lty = 2)


# (c)
x_tilde <- -post$beta[, 2] / (2*post$beta[, 3])

# Plot histogram
hist(x_tilde, col = "#4E79A7",
     main = "",
     xlab = "Time (fraction of year since start)", xlim = c(0.49, 0.52))


# (d) 
# hyperparameters for the 10th-order polynomial regression with shrinkage prior 
p <- 10
mu0_poly <- c(20, -100, 100, rep(0, p - 2))
lambda <- 50
omega0_poly <- lambda * diag(p + 1) 

# Plot prior predictive regression curves for poly
plot(data$time, y, col = "#BAB0AC", pch = 20, main = "", 
     xlab = "Time", ylab = "Temperature", ylim = c(-30, 40))
for (i in 1:N) {
  beta_i <- mvtnorm::rmvnorm(1, mu0_poly, prior_sigma[i] * solve(omega0_poly))
  y_pred <- sapply(time_grid, function(t) sum(beta_i * t^(0:p)))
  lines(time_grid, y_pred, col = "#4E79A7")
}


#------------------------------------------------------------------------------------------------------
# ASSIGNMENT 2

library("mvtnorm")
#(a)
data <- read.csv("Disease.csv")
data <- as.data.frame(data)
#normalize age,symtoms and white blood counts
age_mean <- mean(data$age)
age_sd <- sd(data$age)
duration_mean <- mean(data$duration_of_symptoms)
duration_sd <- sd(data$duration_of_symptoms)
white_mean <- mean(data$white_blood)
white_sd <- sd(data$white_blood)

data[,2] <- (data[,2]-age_mean)/age_sd
data[,4] <- (data[,4]-duration_mean)/duration_sd
data[,6] <- (data[,6]-white_mean)/white_sd

y <- data[,7] #response
X <- as.matrix(data[,1:6]) #Covariates
Xnames <- colnames(X)

#log posterior distribution
LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #log prior
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

#set the prior
mu <- rep(0,ncol(X))
sigma <- 2*diag(ncol(X))
#initial beta
initVal <- rep(0,ncol(X))

#using optim function fine the mode of beta(mean of posterior dist) and observed information at the mode(sigma of posterior dist)
OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,sigma,
                  method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

#find the mu of approx posterior (mode of beta)
approxPostMode <-matrix(OptimRes$par,1,ncol(X))
colnames(approxPostMode) <- Xnames
#find the sd of approx posterior
approxPostStd <- sqrt(diag(solve(-OptimRes$hessian))) # Computing approximate standard deviations.
approxPostStd <- matrix(approxPostStd,1,ncol(X))
colnames(approxPostStd) <- Xnames

#find the 95% equal tail posterior probability interval for coefficient Age
Cred_int <- matrix(0,2,ncol(X)) # Create 95 % approximate credibility intervals for each coefficient
Cred_int[1,] <- approxPostMode - 1.96*approxPostStd
Cred_int[2,] <- approxPostMode + 1.96*approxPostStd
colnames(Cred_int) <- Xnames

cat("95% approximate credibility intervals for Age is [",Cred_int[1,2],Cred_int[2,2],"].","\n")

#Comparison with glm
glmModel<- glm(class_of_diagnosis ~ 0 + ., data = data, family = binomial)
print("Maximum Likelihood Estimates (MLE) for variable Age:")
coef(glmModel)[2]
print("Posterior Means for variable Age:")
approxPostMode[2]

#(b)
#normalizing the sample
xnew <- c(1,38,1,10,0,11000)
xnew_normalized <- c(1,(38-age_mean)/age_sd,1,(10-duration_mean)/duration_sd,
                     0,(11000-white_mean)/white_sd)

#sampling beta from the posterior distribution
set.seed(12345)
beta_sample <- rmvnorm(1000,mean =approxPostMode,sigma=solve(-OptimRes$hessian))

#compute the predictive probability for different beta samples
predict_probs <- numeric(nrow(beta_sample)) 
for (i in 1:nrow(beta_sample)) {
  beta <- beta_sample[i,]
  linPred <- sum(xnew*beta)
  prob <- 1/(1+exp(-linPred))
  predict_probs[i] <- prob
}


hist(predict_probs,breaks = 30,probability = T,main = "Posterior Predictive Distribution",xlab = "P(y=1|x)")