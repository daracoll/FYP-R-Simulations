
```{r}
library(ggplot2)

set.seed(123)

# Parameters
sigma_Z <- 1
sigma_U <- 1
sigma_X <- 1
sigma_Y <- 1
n <- 10000
num_trials <- 100

# Initialize coefficients
true_coefficients <- numeric()
mse_table <- list()
p <- 1

repeat {
  # Generate a new coefficient for the current p
  true_coefficients <- c(true_coefficients, runif(1, min = -5, max = 5))
  
  failed <- FALSE
  mse_list <- matrix(NA, nrow = num_trials, ncol = p)
  estimates <- matrix(NA, nrow = num_trials, ncol = p)
  
  for (trial in 1:num_trials) {
    epsilon_Z <- rnorm(n, mean = 0, sd = sigma_Z)
    U <- rnorm(n, mean = 0, sd = sigma_U)
    epsilon_X <- rnorm(n, mean = 0, sd = sigma_X)
    epsilon_Y <- rnorm(n, mean = 0, sd = sigma_Y)
    
    Z <- epsilon_Z
    X <- Z + U + epsilon_X
    Y <- rowSums(sapply(1:p, function(i) true_coefficients[i] * X^i)) + U + epsilon_Y
    
    A <- matrix(0, nrow = p, ncol = p)
    b <- numeric(p)
    
    for (i in 1:p) {
      for (j in 1:p) {
        A[i, j] <- sum(Z^i * X^j)
      }
      b[i] <- sum(Z^i * Y)
    }
    
    params <- tryCatch(solve(A, b), error = function(e) {
      failed <<- TRUE
      return(rep(NA, p))
    })
    
    if (failed) break
    
    estimates[trial, ] <- params
    mse_list[trial, ] <- (params - true_coefficients[1:p])^2
  }
  
  if (failed) {
    cat("Estimation process stopped for p =", p, "due to singularity.\n")
    break
  }
  
  # Calculate and store MSE for this p
  mse_mean <- colMeans(mse_list, na.rm = TRUE)
  mse_table[[p]] <- mse_mean
  
  # Plot histograms for estimates of each parameter
  for (param_idx in 1:p) {
    param_estimates <- estimates[, param_idx]
    hist_data <- data.frame(Estimate = param_estimates[!is.na(param_estimates)])
    
    p_hist <- ggplot(hist_data, aes(x = Estimate)) +
      geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
      geom_vline(xintercept = true_coefficients[param_idx], color = "red", linetype = "dashed", size = 1) +
      labs(
        title = paste("Histogram of Coefficient", param_idx, "Estimates (p =", p, ")"),
        x = paste("Estimate of Coefficient", param_idx),
        y = "Frequency"
      ) +
      theme_minimal()
    
    print(p_hist)
  }
  
  # Increment p for the next iteration
  p <- p + 1
}

# Format MSE table
mse_df <- data.frame(matrix(NA, nrow = length(mse_table), ncol = max(p - 1, 1)))
colnames(mse_df) <- paste0("Beta", 1:(p - 1))
rownames(mse_df) <- paste0("p=", 1:length(mse_table))

for (i in 1:length(mse_table)) {
  mse_df[i, 1:length(mse_table[[i]])] <- mse_table[[i]]
}

# Display the formatted MSE table
cat("\nFormatted MSE Table:\n")
print(mse_df)

# Optionally save the MSE table as a CSV
 write.csv(mse_df, "formatted_mse_table.csv", row.names = TRUE)

```
