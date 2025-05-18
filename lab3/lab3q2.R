# QUESTION 2

# read data
data <- read.table("D:\\liu\\BL\\Bayesian-Learning\\lab3\\eBayNumberOfBidderData_2025.dat", header = T)

# (a)
mod <- glm(nBids ~ PowerSeller + VerifyID + Sealed + Minblem + MajBlem +
             LargNeg + LogBook + MinBidShare, data = data, family = "poisson")
summary(mod)

# (b)

y <- data$nBids
X <- model.matrix(mod)

# log-posterior function
log_posterior <- function(beta, X, y){
  XtX <- t(X) %*% X
  #sigma0_inv <- 100 * solve(XtX)
  sigma0_inv <- (1/100) * XtX
  # llik 
  eta <- X %*% beta
  llik <- sum(y * eta - exp(eta))
  
  # log-prior (gaussian with mean 0)
  logprior <- -0.5 * t(beta) %*% sigma0_inv %*% beta
  
  # Return scalar value
  log_post <- llik + logprior
    
  return(as.numeric(log_post))
}

# find posterior mode
# Initialize beta
init_beta <- rep(0, ncol(X))

# Optimization
opt_result <- optim(
  par = init_beta,
  fn = function(beta) -log_posterior(beta, X, y), # negative because optim minimizes
  method = "BFGS",
  hessian = TRUE
)

# Posterior mode and Hessian
beta_tilde <- opt_result$par
Hessian <- opt_result$hessian


# (c)

MetropolisHastings <- function(x_init, n_iter, prop_dist_rand, target, X, y, Sigma, c){

  # initial values
  p <- length(x_init)
  x <- matrix(NA, nrow = n_iter, ncol = p)
  x[1, ] <- x_init
  acc <- 0
  
  for (t in 2:n_iter) {
    x_candidate <- prop_rand(x[t - 1, ], Sigma, c)
    x_candidate <- as.vector(x_candidate)
    
    log_r <- target(x_candidate, X, y) - target(x[t - 1, ], X, y)
    if (log(runif(1)) < log_r) {
      x[t, ] <- x_candidate
      acc <- acc + 1
    } else {
      x[t, ] <- x[t - 1, ]
    }
  }
  return(list(chain = x, acceptance_rate = sum(acc)/(n_iter-1)))
}

prop_rand <- function(theta, Sigma, c){
  mvtnorm::rmvnorm(1, mean = theta, sigma = c * Sigma)
}

#init_beta <- rep(0, ncol(X))  # or use beta_tilde from earlier optimization
Sigma <- solve(Hessian)

result <- MetropolisHastings(
  x_init = beta_tilde, # init_beta,
  n_iter = 10000,
  prop_dist = prop_rand,
  target = log_posterior,
  X = X,
  y = y,
  Sigma = Sigma,
  c = 0.5  # tuning parameter
)


posterior_samples <- result$chain
acceptance_rate <- result$acceptance_rate

print(acceptance_rate)
plot.ts(posterior_samples[,1], main = "Traceplot for β₀")
hist(posterior_samples[,1], main = "Posterior of β₀", breaks = 40)

par(mfrow = c(3, 3))
for (j in 1:ncol(posterior_samples)) {
  plot(posterior_samples[, j], type = "l", main = paste("Traceplot: Beta", j),
       xlab = "Iteration", ylab = paste("β[", j, "]", sep = ""))
}

# Density plots
par(mfrow = c(3, 3))
for (j in 1:ncol(posterior_samples)) {
  plot(density(posterior_samples[, j]), main = paste("Density: Beta", j))
}

# (d)

x_new <- c(1,    # Intercept
           1,    # PowerSeller
           0,    # VerifyID
           1,    # Sealed
           0,    # MinBlem
           1,    # MajBlem
           0,    # LargNeg
           1.3,  # LogBook
           0.7)  # MinBidShare

# Compute lambda for each posterior draw
lambda_draws <- apply(posterior_samples, 1, function(beta) {
  exp(sum(beta * x_new))
})

# Simulate nBids for new auction
set.seed(1234)
pred_nBids <- rpois(length(lambda_draws), lambda = lambda_draws)

# Plot predictive distribution
par(mfrow = c(1,1))
hist(pred_nBids, col = "skyblue",
     main = "Predictive Distribution of Number of Bidders",
     xlab = "Number of Bidders", freq = FALSE)

# Estimate P(nBids = 0)
prob_zero_bidders <- mean(pred_nBids == 0)
cat("Estimated probability of zero bidders:", prob_zero_bidders, "\n")

