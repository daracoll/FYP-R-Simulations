```{r}
library(ggplot2)

set.seed(348)

alpha <- 2
sample_sizes <- c(10, 100, 1000, 10000,100000)
n_reps <- 10000
results <- matrix(NA, nrow = length(sample_sizes), ncol = 2)

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

    A <- sum(Z * X)
    b <- sum(Z * Y)

    alpha_hat <- b / A
    alpha_estimates[j] <- alpha_hat
  }
  
  mse_alpha <- mean((alpha_estimates - alpha)^2)
  var_alpha <- var(alpha_estimates)

  results[i, 1] <- mse_alpha
  results[i, 2] <- var_alpha

  hist_data_alpha <- data.frame(Estimate = alpha_estimates)

  p_hist_alpha <- ggplot(hist_data_alpha, aes(x = Estimate)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = alpha, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Histogram of Alpha Estimates (n =", n, ")"),
      x = "Estimate of Alpha",
      y = "Frequency"
    ) +
    theme_minimal()

  print(p_hist_alpha)
}

colnames(results) <- c("MSE", "Variance")
rownames(results) <- as.character(sample_sizes)
print(results)

```
