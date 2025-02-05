---
title: "what i need"
output: html_document
date: "2025-02-03"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```



```{r}
library(ggplot2)
library(AER)
library(dHSIC)
library(dplyr)
library(doParallel)
library(foreach)

set.seed(42)

# Simulation parameters
n <- 500  # Sample size
num_trials <- 100  # Number of trials
true_beta <- 11
sigma_U <- 2
sigma_X <- 2
sigma_Y <- 3

# Function to generate Gaussian Z
simulate_Z <- function(n) {
  return(rnorm(n))
}

# Function to estimate using OLS estimator (assumes no intercept)
ols_estimator <- function(X, Y) {
  theta_ols <- tryCatch(coef(lm(Y ~ X - 1))[1], error = function(e) NA)
  if (is.na(theta_ols)) { theta_ols <- 1 }
  return(theta_ols)
}

# Function to estimate using 2SLS IV estimator (assumes no intercept)
iv_estimator <- function(Z, X, Y) {
  theta_iv <- tryCatch(coef(ivreg(Y ~ X - 1 | Z))[1], error = function(e) NA)  
  if (is.na(theta_iv)) { theta_iv <- ols_estimator(X, Y) }
  return(theta_iv)
}

# HSIC-X Gradient Function
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

# HSIC-X Estimator with Restarting and Convergence Check
hsicx_estimator <- function(X, Y, Z, kernel_type, t = 20, k = 100, m = 32, gamma = 0.1, alpha = .05) {
  l <- 0
  pval <- 0
  update <- 1
  theta_iv <- iv_estimator(Z, X, Y)  # Store IV estimate as fallback
  theta_best <- theta_iv  # Default to IV estimate
  pval_best <- 0  # Track best p-value
  
  while (l < t && pval < alpha) {
    if (l == 0) {
      theta <- theta_iv  # Initialize with IV estimate
    } else {
      theta <- runif(1, min = 7, max = 15)  # Random restart
    }
    
    while (update > 0.001) {
      theta_old <- theta
      res <- Y - (theta * X)
      
      for (i in 1:k) {
        sample <- sample_n(data.frame(X, Y, Z), m)
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
    
    # Update best theta if new p-value is higher
    if (pval > pval_best) {
      pval_best <- pval
      theta_best <- theta
    }
    
    # Print progress update
    cat("Iteration:", l, "Theta Estimate:", theta, "P-value:", pval, "\n")
  }
  
  # If t is exceeded, return best theta found instead of IV estimate
  return(theta_best)
}

# Run simulation
run_simulation <- function(num_trials, n) {
  estimates_list <- list()
  
  for (trial in 1:num_trials) {
    # Generate instruments
    Z_cont <- simulate_Z(n)
    
    # Generate error terms
    U <- rnorm(n, 0, sigma_U)
    epsilon_X <- rnorm(n, 0, sigma_X)
    
    # Continuous Instrument Data
    X_cont <- Z_cont + U + epsilon_X  
    Y_cont <- true_beta * X_cont + U + rnorm(n, 0, sigma_Y)
    X_nl_cont <- Z_cont * U -2* U + epsilon_X
    Y_nl_cont <- true_beta * X_nl_cont +5* U + rnorm(n, 0, sigma_Y)
    
    # OLS Estimates
    beta_ols_cont <- ols_estimator(X_cont, Y_cont)
    beta_ols_nl_cont <- ols_estimator(X_nl_cont, Y_nl_cont)
    
    # IV Estimates
    beta_iv_cont <- iv_estimator(Z_cont, X_cont, Y_cont)
    beta_iv_nl_cont <- iv_estimator(Z_cont, X_nl_cont, Y_nl_cont)
    
    # HSIC-X Estimates
    beta_hsicx_cont <- hsicx_estimator(X_cont, Y_cont, Z_cont, "gaussian")
    beta_hsicx_nl_cont <- hsicx_estimator(X_nl_cont, Y_nl_cont, Z_cont, "gaussian")
    
    # Store estimates
    estimates_list[[trial]] <- list(
      Trial = trial,
      beta_ols_cont = beta_ols_cont, beta_iv_cont = beta_iv_cont, beta_hsicx_cont = beta_hsicx_cont,
      beta_ols_nl_cont = beta_ols_nl_cont, beta_iv_nl_cont = beta_iv_nl_cont, beta_hsicx_nl_cont = beta_hsicx_nl_cont
    )
  }
  return(estimates_list)
}

# Run and process results
results <- run_simulation(num_trials, n)
results_estimation_data <- bind_rows(lapply(results, as.data.frame))

# Print stored estimates
print(results_estimation_data)
write.csv(results_estimation_data, "HSIC-X Estimation Data", row.names = FALSE)

```


```{r}
library(tidyr)
true_value <- 11

# Calculate MSE for each method (including OLS)
mse_data <- results_estimation_data %>%
  summarise(
    mse_ols_cont = mean((beta_ols_cont - true_value)^2),
    mse_iv_cont = mean((beta_iv_cont - true_value)^2),
    mse_hsicx_cont = mean((beta_hsicx_cont - true_value)^2),
    mse_ols_nl_cont = mean((beta_ols_nl_cont - true_value)^2),
    mse_iv_nl_cont = mean((beta_iv_nl_cont - true_value)^2),
    mse_hsicx_nl_cont = mean((beta_hsicx_nl_cont - true_value)^2)
  ) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "MSE")

# Separate MSE data for linear and non-linear cases
mse_data <- mse_data %>%
  mutate(Type = ifelse(grepl("nl", Method), "Non-Linear", "Linear"))

# Plot Average MSE for Linear case
mse_plot_linear <- ggplot(mse_data %>% filter(Type == "Linear"), aes(x = Method, y = MSE, fill = Method)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Average MSE for Linear Case", x = "Method", y = "Mean Squared Error")

print(mse_plot_linear)

# Plot Average MSE for Non-Linear case
mse_plot_nonlinear <- ggplot(mse_data %>% filter(Type == "Non-Linear"), aes(x = Method, y = MSE, fill = Method)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Average MSE for Non-Linear Case", x = "Method", y = "Mean Squared Error")

print(mse_plot_nonlinear)

# Function to plot histograms with true value as reference
hist_plot <- function(column_name, title) {
  ggplot(results_estimation_data, aes_string(x = column_name)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = true_value, linetype = "dashed", color = "red", size = 1) +
    theme_minimal() +
    labs(title = title, x = "Estimate", y = "Density")
}

# Generate histograms for each method (including OLS)
hist_ols_cont <- hist_plot("beta_ols_cont", "Histogram of OLS Estimates (Continuous)")
hist_iv_cont <- hist_plot("beta_iv_cont", "Histogram of IV Estimates (Continuous)")
hist_hsicx_cont <- hist_plot("beta_hsicx_cont", "Histogram of HSIC-X Estimates (Continuous)")
hist_ols_nl_cont <- hist_plot("beta_ols_nl_cont", "Histogram of OLS Estimates (Non-Linear Continuous)")
hist_iv_nl_cont <- hist_plot("beta_iv_nl_cont", "Histogram of IV Estimates (Non-Linear Continuous)")
hist_hsicx_nl_cont <- hist_plot("beta_hsicx_nl_cont", "Histogram of HSIC-X Estimates (Non-Linear Continuous)")

# Print histograms
print(hist_ols_cont)
print(hist_iv_cont)
print(hist_hsicx_cont)
print(hist_ols_nl_cont)
print(hist_iv_nl_cont)
print(hist_hsicx_nl_cont)

```


