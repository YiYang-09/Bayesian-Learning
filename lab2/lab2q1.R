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

scale_inv_chi <- function(nu0, sigma0, N) {
   sigma0 <- as.numeric(sigma0)
   samples <- (nu0 * sigma0) / rchisq(N, df = nu0)
   return(samples)
 }

N <- 50
set.seed(123)
prior_sigma <- scale_inv_chi(nu0, sigma0, N)

prior_beta <- matrix(NA, nrow = N, ncol = 3)
for(i in 1:N){
  prior_beta[i, ] <- mvtnorm::rmvnorm(n = 1, mean = mu0, sigma = prior_sigma[i] * solve(omega0))
}

# Plot prior regression curves
plot(data$time, y, col = "#BAB0AC", pch = 20, main = "Prior Predictive Regression Curves", 
     xlab = "Time", ylab = "Temperature", ylim = c(-30, 40))
# plot(1, type = "n", main = "Prior Predictive Regression Curves", 
#      xlab = "Time", ylab = "Temperature", ylim = c(-30, 40), xlim = c(0, 1))
for (i in 1:N) {
  curve <- prior_beta[i, 1] + prior_beta[i, 2] * data$time + prior_beta[i, 3] * data$time^2
  lines(data$time, curve, col = "#4E79A7")
}


# (b)
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
hist(post$sigma, col = "#4E79A7", 
     main = expression(paste("Marginal posterior of ", sigma)), xlab = expression(sigma))
hist(post$beta[, 1], col = "#4E79A7",
     main = expression(paste("Marginal posterior of ", beta[0])), xlab = expression(beta[0]))
hist(post$beta[, 2], col = "#4E79A7", 
     main = expression(paste("Marginal posterior of ", beta[1])), xlab = expression(beta[1]))
hist(post$beta[, 3], col = "#4E79A7",
     main = expression(paste("Marginal posterior of ", beta[2])), xlab = expression(beta[2]))

time_grid <- seq(min(data$time), max(data$time), length.out = 100)
f_draws <- matrix(NA, nrow = N, ncol = length(time_grid))
for (i in 1:N) {
  f_draws[i, ] <- post$beta[i, 1] + post$beta[i, 2] * time_grid + post$beta[i, 3] * time_grid^2
}

draws_median <- apply(f_draws, 2, median)
draws_lb <- apply(f_draws, 2, quantile, probs = 0.05)
draws_ub <- apply(f_draws, 2, quantile, probs = 0.95)

# Plot data
plot(x, y, col = "#BAB0AC", pch = 20,
     main = "Posterior Regression Curve with 90% Credible Intervals",
     xlab = "Time", ylab = "Temperature")

# Add posterior median curve
lines(time_grid, draws_median, col = "#4E79A7", lwd = 2)

# Add 90% posterior credible interval
lines(time_grid, draws_lb, col = "#E15759", lty = 2)
lines(time_grid, draws_ub, col = "#E15759", lty = 2)


# (c)

x_tilde <- -post$beta[, 2] / (2*post$beta[, 3])

# Plot histogram
hist(x_tilde, col = "#4E79A7",
     main = "Posterior Distribution of x~ (Minimum Temperature Time)",
     xlab = "Time (fraction of year since start)", xlim = c(0.49, 0.52))


# (d) 10th-order polynomial regression with shrinkage prior
p <- 10
mu0_poly <- c(20, -100, 100, rep(0, p - 2))
lambda <- 50
omega0_poly <- lambda * diag(p + 1) 

# Plot prior predictive regression curves for poly
plot(data$time, y, col = "#BAB0AC", pch = 20, main = "Prior Predictive Regression Curves", 
     xlab = "Time", ylab = "Temperature", ylim = c(-30, 40))
for (i in 1:N) {
  beta_i <- mvtnorm::rmvnorm(1, mu0_poly, prior_sigma[i] * solve(omega0_poly))
  y_pred <- sapply(time_grid, function(t) sum(beta_i * t^(0:p)))
  lines(time_grid, y_pred, col = "#4E79A7")
}

