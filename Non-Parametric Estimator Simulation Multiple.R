


```{r}
library(ggplot2)

set.seed(123)

sigma_Z <- 1
sigma_U <- 1
sigma_X <- 1
sigma_Y <- 1
n <- 10000
num_trials <- 100

true_coefficients <- numeric()
mse_table <- list()
p <- 1

repeat {
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
  
  if (failed) break
  
  mse_mean <- colMeans(mse_list, na.rm = TRUE)
  mse_table[[p]] <- mse_mean
  
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
  
  p <- p + 1
}

mse_df <- data.frame(matrix(NA, nrow = length(mse_table), ncol = max(p - 1, 1)))
colnames(mse_df) <- paste0("Beta", 1:(p - 1))
rownames(mse_df) <- paste0("p=", 1:length(mse_table))

for (i in 1:length(mse_table)) {
  mse_df[i, 1:length(mse_table[[i]])] <- mse_table[[i]]
}

cat("\nFormatted MSE Table:\n")
print(mse_df)


```




```{r}
library(ggplot2)
library(tidyr)
library(dplyr)

mse_df <- data.frame(matrix(NA, nrow = length(mse_table), ncol = max(p - 1, 1)))
colnames(mse_df) <- paste0("Beta", 1:(p - 1))
rownames(mse_df) <- paste0("p=", 1:length(mse_table))

for (i in 1:length(mse_table)) {
  mse_df[i, 1:length(mse_table[[i]])] <- mse_table[[i]]
}

mse_long <- mse_df %>%
  rownames_to_column("p") %>%
  gather(key = "Parameter", value = "Observed_MSE", -p)

mse_long$p <- as.integer(gsub("p=", "", mse_long$p))

parameter_to_plot <- "Beta2"

filtered_mse_results <- mse_long[mse_long$Parameter == parameter_to_plot, ]

mse_plot <- ggplot(filtered_mse_results, aes(x = p, y = Observed_MSE)) +
  geom_line(size = 1, colour = "blue") +
  geom_point(size = 2, colour = "blue") +
  labs(
    title = paste("MSE of", parameter_to_plot, "vs. p"),
    x = "p (Number of Parameters)",
    y = "Observed MSE"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = filtered_mse_results$p) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(mse_plot)




```

 






