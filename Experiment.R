---
title: "IV Experiment"
output: html_document
date: "2025-02-04"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```









```{r}
# Load required packages

library(AER)      
library(dHSIC)   
library(dplyr)    
library(lmtest)   
library(sandwich) 
library(kernlab)  

# Load the CigarettesSW dataset
data("CigarettesSW", package = "AER")  
CigarettesSW <- as.data.frame(na.omit(CigarettesSW))  # Remove NAs

# Define variables
Y_cig <- CigarettesSW$packs  # Dependent variable: Cigarette consumption per capita
X_cig <- CigarettesSW$price  # Endogenous variable: Price of cigarettes
Z_cig <- CigarettesSW$tax    # Instrument: State-level excise tax

# --- 🔹 First-Stage Regression (Checking Instrument Strength) ---
first_stage_cig <- lm(X_cig ~ Z_cig, data = CigarettesSW)
summary(first_stage_cig)

# Extract F-statistic for weak instrument check
first_stage_fstat_cig <- summary(first_stage_cig)$fstatistic[1]
cat("First-Stage F-Statistic (CIGARETTES):", first_stage_fstat_cig, "\n")

# --- 🔹 2SLS Estimation (Accounting for Intercept) ---
iv_model_cig <- ivreg(Y_cig ~ X_cig | Z_cig, data = CigarettesSW)
summary(iv_model_cig)
theta_2sls_cig <- coef(iv_model_cig)  # Includes both intercept and coefficient

# --- 🔹 HSIC-X Gradient Function (Now With Intercept) ---
hsic_grad <- function(x, y, z, theta){
  theta1 <- theta[2] - 0.001
  theta2 <- theta[2] + 0.001
  resid1 <- y - (theta[1] + theta1 * x)
  resid2 <- y - (theta[1] + theta2 * x)
  hsic1 <- dhsic(matrix(resid1, ncol = 1), matrix(z, ncol = 1))$dHSIC
  hsic2 <- dhsic(matrix(resid2, ncol = 1), matrix(z, ncol = 1))$dHSIC
  return((hsic2 - hsic1) / (theta2 - theta1))
}

# HSIC-X Estimator with Fixed Intercept and Early Stopping
hsicx <- function(X, Y, Z, t = 100, k = 100, m = 32, gamma = 0.1, alpha = .05) {
  theta_best <- c(theta_2sls_cig[1], theta_2sls_cig[2])  # Initial theta using IV estimates
  pval_best <- 0
  l <- 0  # Initialize iteration counter
  pval <- 0  # Initialize p-value
  
  while (l < t && pval < alpha) {
    if (l == 0) {
      theta <- c(theta_2sls_cig[1], theta_2sls_cig[2])  # Initialize with 2SLS
    } else {
      # Sample intercept and slope from normal distributions centered at IV estimates
      theta <- c(rnorm(1, mean = theta_2sls_cig[1], sd = 1), 
                 rnorm(1, mean = theta_2sls_cig[2], sd = 1))
    }
    
    update <- 1  # Reset update value
    while (update > 0.001) {  # Ensure updates occur until convergence
      theta_old <- theta[2]
      
      for (i in 1:k) {
        sample_data <- slice_sample(data.frame(X, Y, Z), n = min(m, nrow(CigarettesSW)))  # Random subsample
        theta[2] <- theta[2] - gamma * hsic_grad(sample_data$X, sample_data$Y, sample_data$Z, theta)
      }
      
      update <- sqrt((theta_old - theta[2])^2)  # Compute update size
    }
    
    resid <- matrix(Y - (theta[1] + theta[2] * X), ncol = 1)  # Adjusted for new intercept
    test <- dhsic.test(resid, matrix(Z, ncol = 1))  # Test dependence
    pval <- test$p.value  # Extract p-value
    
    # **Print Iteration Updates**
    cat("Iteration:", l + 1, "| Intercept:", theta[1], "| Coefficient:", theta[2], "| P-value:", pval, "\n")
    
    # If p-value exceeds alpha, stop and return the current theta
    if (pval > alpha) {
      cat("Stopping early: P-value exceeded alpha at iteration", l + 1, "\n")
      return(theta)  # Exit function immediately
    }
    
    if (pval > pval_best) {  # Update best values
      pval_best <- pval
      theta_best <- theta
    }
    
    l <- l + 1  # Increment iteration counter
  }
  
  return(theta_best)
}


# --- 🔹 Run HSIC-X Estimation ---
theta_hsicx_cig <- hsicx(X_cig, Y_cig, Z_cig)

# --- 🔹 Print Results ---
cat("2SLS Estimate (CIGARETTES): Intercept:", theta_2sls_cig[1], "Coefficient:", theta_2sls_cig[2], "\n")
cat("HSIC-X Estimate (CIGARETTES): Intercept:", theta_hsicx_cig[1], "Coefficient:", theta_hsicx_cig[2], "\n")

```




```{r}
# --- 🔹 Bootstrapping for 2SLS Confidence Interval ---
set.seed(123)  # For reproducibility
n_boot <- 100
boot_2sls <- numeric(n_boot)

for (i in 1:n_boot) {
  sample_indices <- sample(1:nrow(CigarettesSW), replace = TRUE)
  sample_data <- CigarettesSW[sample_indices, ]
  
  # Re-estimate IV model on bootstrap sample
  boot_model <- ivreg(packs ~ price | tax, data = sample_data)
  boot_2sls[i] <- coef(boot_model)[2]  # Store coefficient estimate
}

# Compute 95% Confidence Interval for 2SLS
ci_2sls <- quantile(boot_2sls, probs = c(0.025, 0.975))

# --- 🔹 Bootstrapping for HSIC-X Confidence Interval ---
set.seed(123)
boot_hsicx <- numeric(n_boot)

for (i in 1:n_boot) {
  sample_indices <- sample(1:nrow(CigarettesSW), replace = TRUE)
  sample_data <- CigarettesSW[sample_indices, ]
  
  # Run HSIC-X estimation on bootstrap sample
  boot_hsicx[i] <- hsicx(sample_data$price, sample_data$packs, sample_data$tax)[2]  # Store estimated coefficient
}

# Compute 95% Confidence Interval for HSIC-X
ci_hsicx <- quantile(boot_hsicx, probs = c(0.025, 0.975))

# --- 🔹 Print Results ---
cat("2SLS Estimate:", theta_2sls_cig[2], "\n")
cat("2SLS 95% CI:", ci_2sls, "\n")
cat("HSIC-X Estimate:", theta_hsicx_cig[2], "\n")
cat("HSIC-X 95% CI:", ci_hsicx, "\n")

```

