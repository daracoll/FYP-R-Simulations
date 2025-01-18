
```{r}
library(ggplot2)

set.seed(348)

true_beta_1 <- 1
true_beta_2 <- 4
sigma_Z <- 1
sigma_U <- 1
sigma_X <- 1
sigma_Y <- 1
sample_sizes <- c(10, 100, 1000, 10000)
n_reps <- 100
mse_results <- matrix(NA, nrow = length(sample_sizes), ncol = 2)

for (i in 1:length(sample_sizes)) {
  
  n <- sample_sizes[i]
  beta_1_estimates <- numeric(n_reps)
  beta_2_estimates <- numeric(n_reps)
  
  for (j in 1:n_reps) {
    epsilon_Z <- rnorm(n, mean = 0, sd = sigma_Z)
    U <- rnorm(n, mean = 0, sd = sigma_U)
    epsilon_X <- rnorm(n, mean = 0, sd = sigma_X)
    epsilon_Y <- rnorm(n, mean = 0, sd = sigma_Y)

    Z <- epsilon_Z
    X <- Z + U + epsilon_X
    Y <- true_beta_1 * X + true_beta_2 * X^2 + U + epsilon_Y

    data <- data.frame(Z = Z, X = X, Y = Y)

    A11 <- sum(Z * X)
    A12 <- sum(Z * X^2)
    A21 <- sum(Z^2 * X)
    A22 <- sum(Z^2 * X^2)

    b1 <- sum(Z * Y)
    b2 <- sum(Z^2 * Y)

    A <- matrix(c(A11, A12, A21, A22), nrow = 2, byrow = TRUE)
    b <- c(b1, b2)

    params <- solve(A, b)
    beta_1_estimates[j] <- params[1]
    beta_2_estimates[j] <- params[2]
  }
  
  mse_beta_1 <- mean((beta_1_estimates - true_beta_1)^2)
  mse_beta_2 <- mean((beta_2_estimates - true_beta_2)^2)

  mse_results[i, 1] <- mse_beta_1
  mse_results[i, 2] <- mse_beta_2

  hist_data_beta_1 <- data.frame(Estimate = beta_1_estimates)
  hist_data_beta_2 <- data.frame(Estimate = beta_2_estimates)

  p_hist_beta_1 <- ggplot(hist_data_beta_1, aes(x = Estimate)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = true_beta_1, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Histogram of beta_1 Estimates (n =", n, ")"),
      x = "Estimate of beta_1",
      y = "Frequency"
    ) +
    theme_minimal()

  p_hist_beta_2 <- ggplot(hist_data_beta_2, aes(x = Estimate)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = true_beta_2, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Histogram of beta_2 Estimates (n =", n, ")"),
      x = "Estimate of beta_2",
      y = "Frequency"
    ) +
    theme_minimal()

  print(p_hist_beta_1)
  print(p_hist_beta_2)
}

colnames(mse_results) <- c("beta_1", "beta_2")
rownames(mse_results) <- as.character(sample_sizes)
print(mse_results)

```
