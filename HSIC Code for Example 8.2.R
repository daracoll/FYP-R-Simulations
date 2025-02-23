```{r}
library(ggplot2)
library(dHSIC)
library(dplyr)

set.seed(348)

alpha <- 2
sample_sizes <- c(10,100,1000,10000)
n_reps <- 50
results <- matrix(NA, nrow = length(sample_sizes), ncol = 2)

all_estimates <- list()

hsic_grad <- function(x, y, z, theta, kernel_type) {
  theta1 <- theta - 0.001
  theta2 <- theta + 0.001
  resid1 <- y - (theta1 * x)
  resid2 <- y - (theta2 * x)
  hsic1 <- dhsic(resid1, z, kernel = kernel_type)$dHSIC
  hsic2 <- dhsic(resid2, z, kernel = kernel_type)$dHSIC
  slope <- (hsic2 - hsic1) / (theta2 - theta1)
  return(slope)
}

hsicx_estimator <- function(X, Y, Z, kernel_type, t = 30, k = 100, m = 32, gamma = 0.1, alpha = .05) {
  l <- 0
  pval <- 0
  update <- 1
  theta <- runif(1, min = 0, max = 4)
  theta_best <- theta
  pval_best <- 0  
  
  while (l < t && pval < alpha) {
    while (update > 0.001) {
      theta_old <- theta
      res <- Y - (theta * X)
      
      for (i in 1:k) {
        sample_size <- min(m, length(X))
        sample <- sample_n(data.frame(X, Y, Z), sample_size, replace = (sample_size < m))
        res_sample <- sample$Y - (theta * sample$X)
        loss <- dhsic(res_sample, sample$Z, kernel = kernel_type)
        theta <- theta - gamma * hsic_grad(sample$X, sample$Y, sample$Z, theta, kernel_type)
      }
      
      update <- sqrt((theta_old - theta)^2)
    }
    
    resid <- Y - (theta * X)
    test <- dhsic.test(resid, Z, kernel = kernel_type)
    pval <- test$p.value
    l <- l + 1
    
    if (pval > pval_best) {
      pval_best <- pval
      theta_best <- theta
    }
  }
  
  return(theta_best)
}

for (i in 1:length(sample_sizes)) {
  n <- sample_sizes[i]
  alpha_estimates <- numeric(n_reps)
  
  for (j in 1:n_reps) {
    epsilon_Z <- rnorm(n)
    epsilon_U <- rnorm(n)
    epsilon_X <- rnorm(n)
    epsilon_Y <- rnorm(n)

    Z <- epsilon_Z
    U <- epsilon_U
    X <- Z * U + epsilon_X
    Y <- alpha * X + U + epsilon_Y

    alpha_hat <- hsicx_estimator(X, Y, Z, "gaussian")
    alpha_estimates[j] <- alpha_hat
  }
  
  mse_alpha <- mean((alpha_estimates - alpha)^2)
  var_alpha <- var(alpha_estimates)

  results[i, 1] <- mse_alpha
  results[i, 2] <- var_alpha

  all_estimates[[i]] <- data.frame(Estimate = alpha_estimates, Sample_Size = n)
}

all_estimates_df <- bind_rows(all_estimates)

colnames(results) <- c("MSE", "Variance")
rownames(results) <- as.character(sample_sizes)

print(results)
print(all_estimates_df)

write.csv(all_estimates_df, "hsic_estimates.csv", row.names = FALSE)
```
