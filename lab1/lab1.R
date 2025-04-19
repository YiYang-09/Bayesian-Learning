# ASSIGNMENT 1
# data and parameters
n <- 78
f <- 35
s <- n - f
alpha0 <- 7
beta0 <- 7
alpha <- alpha0 + s
beta <- beta0 + f

# calculate theoretical mean and sd for the beta posterior
true_mean <- alpha / (alpha + beta)
true_sd <- sqrt(alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1)))

# (a)
# sample of size 10000 from the beta posterior distribution
set.seed(12345)
draw <- rbeta(10000, shape1 = alpha, shape2 = beta)

# calculate mean and sd for each draw
mean_vals <- numeric(9999)
sd_vals <- numeric(9999)
for(i in 1:10000){
  mean_vals[i-1] <- mean(draw[1:i])
  sd_vals[i-1] <- sd(draw[1:i])
}

# mean values plot with theoretical value
ggplot(mapping = aes(x = seq_along(mean_vals), y = mean_vals)) +
  geom_line(color = "#4E79A7") +
  geom_hline(yintercept = true_mean, color = "#E15759", linetype = "dashed") +
  labs(title = "Convergence of Posterior Mean",
       x = "Sample Size",
       y = "Estimated Mean") +
  theme_minimal()

# sd values plot with theoretical value
ggplot(mapping = aes(x = seq_along(sd_vals), y = sd_vals)) +
  geom_line(color = "#4E79A7") +
  geom_hline(yintercept = true_sd, color = "#E15759", linetype = "dashed") +
  labs(title = "Convergence of Posterior Standard Deviation",
       x = "Sample Size",
       y = "Estimated SD") +
  theme_minimal()

# (b)
# sample probability
prob <- sum((draw > 0.5))/10000 

# theoretical probability with posterior distribution
true_prob <- pbeta(0.5, alpha, beta, lower.tail = F)

# (c)
# calculate odds
odds <- draw/(1 - draw)

# plot posterior distribution of odds
ggplot(mapping = aes(x = odds)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "#4E79A7", color = "white", alpha = 0.7) +
  geom_density(color = "#E15759", size = 1) +
  labs(title = expression(paste("Distribution of ", phi)),
       x = expression(phi),
       y = "Density") +
  theme_minimal()


# ASSIGNMENT 2
# data and parameters
y <- c(22, 33, 31, 29, 65, 78, 17, 24)
n <- 8
mu <- 3.65
tau <- sum((log(y) - mu)^2)/n

# (a)
# posterior distribution of sigma^2 (scale-inv-chi)
set.seed(12345)
samples <- (n * tau) / rchisq(10000, df = n)

# plot posterior distribution of sigma^2
ggplot(mapping = aes(x = samples)) +
  geom_histogram(bins = 50, fill = "#4E79A7", color = "white", alpha = 0.8) +
  labs(
    title = expression(paste("Posterior of ", sigma^2)),
    x = expression(sigma^2),
    y = "Frequency"
  ) +
  theme_minimal()

# (b)
# calculate gini index with the provided formula
gini <- 2 * pnorm(sqrt(samples)/sqrt(2), 0, 1) - 1

# plot posterior distribution of gini index
ggplot(mapping = aes(x = gini)) +
  geom_histogram(bins = 50, fill = "#4E79A7", color = "white", alpha = 0.8) +
  labs(
    title = "Posterior of G",
    x = expression(G),
    y = "Frequency"
  ) +
  theme_minimal()

# (c)
# compute equal-tail credible interval
ci <- quantile(gini, c(0.025, 0.975))

# (d)
# compute HPDI
hpdi <- hdi(gini, ci = 0.95)

table <- data.frame("Lower bound" = c(ci[1], hpdi[1,2]), "Upper bound" = c(ci[2], hpdi[1,3]))
rownames(table) <- c("Equal tail", "HPDI")
kable(table, caption = "Confidence intervals for G")


# ASSIGNMENT 3
# (a)
# data and parameters
y <- c(0, 2, 5, 5, 7, 1, 4)
sigma <- 5

# the expression of posterior distribution unnormalized
posterior <- function(lambda, sigma, y){
  s <- sum(y)
  n <- length(y)
  res <- lambda^(s)*exp(-n*lambda - (lambda^2)/(2*sigma^2))
  return(res)
}

# unnormalized posterior distribution
lambda_grid <- seq(0, 7, by=0.05)
unnormalized <- posterior(lambda_grid, sigma, y)

# normalized posterior
normalized <- unnormalized / sum(unnormalized * 0.05)

# plot normalized posterior
df <- data.frame(lambda = lambda_grid, density = normalized)
ggplot(data = df, mapping = aes(x = lambda, y = density))+
  geom_point(color = "#4E79A7")+
  labs(title = "Plot of posterior distribution ")+
  theme_classic()

#(b)
mode_index <- which.max(normalized)
mode <- lambda_grid[mode_index]
