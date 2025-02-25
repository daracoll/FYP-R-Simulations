```{r}
library(ggplot2)
library(AER)
library(dHSIC)
library(dplyr)

set.seed(42)

n <- 500  
num_trials <- 50  
true_beta <- 11
sigma_U <- 3
sigma_X <- 4
sigma_Y <- 2


simulate_Z_cont <- function(n) {
  return(rnorm(n))
}

simulate_Z_disc <- function(n) {
  return(sample(1:3, n, replace = TRUE))  
}


ols_estimator <- function(X, Y) {
  theta_ols <- tryCatch(coef(lm(Y ~ X - 1))[1], error = function(e) NA)
  if (is.na(theta_ols)) { theta_ols <- 1 }
  return(theta_ols)
}


iv_estimator <- function(Z, X, Y) {
  Z_numeric <- as.numeric(Z)
  theta_iv <- tryCatch(coef(ivreg(Y ~ X - 1 | Z_numeric))[1], error = function(e) NA)  
  if (is.na(theta_iv)) { theta_iv <- ols_estimator(X, Y) }
  return(theta_iv)
}


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
  theta_iv <- iv_estimator(Z, X, Y)  
  theta_best <- theta_iv  
  pval_best <- 0  
 
  while (l < t && pval < alpha) {
    if (l == 0) {
      theta <- theta_iv  
    } else {
      theta <- runif(1, min = 8, max = 14)  
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
   
    if (pval > pval_best) {
      pval_best <- pval
      theta_best <- theta
    }
   
    cat("Iteration:", l, "Theta Estimate:", theta, "P-value:", pval, "\n")
  }
 
  if (l == t && pval < alpha) {
    warning(" HSIC-X estimator did not converge within ", t,
            " iterations. Using best estimate found: ", round(theta_best, 4))
  }
 
  return(theta_best)
}



run_simulation <- function(num_trials, n) {
  results <- list()
 
  for (trial in 1:num_trials) {
    Z_cont <- simulate_Z_cont(n)
    Z_disc <- simulate_Z_disc(n)
   
    U <- rnorm(n, 0, sigma_U)
    epsilon_X <- rnorm(n, 0, sigma_X)
   
    X_cont <- Z_cont -2* U + epsilon_X  
    Y_cont <- true_beta * X_cont +5* U + rnorm(n, 0, sigma_Y)
    X_nl_cont <- Z_cont * U -2*U + epsilon_X
    Y_nl_cont <- true_beta * X_nl_cont +  5*U + rnorm(n, 0, sigma_Y)
   
    X_disc <- Z_disc -2* U + epsilon_X  
    Y_disc <- true_beta * X_disc +5* U + rnorm(n, 0, sigma_Y)
    X_nl_disc <-  (Z_disc ) * U  -2*U+ epsilon_X
    Y_nl_disc <- true_beta * X_nl_disc + 5* U + rnorm(n, 0, sigma_Y)
   
   
    beta_ols_cont <- ols_estimator(X_cont, Y_cont)
    beta_iv_cont <- iv_estimator(Z_cont, X_cont, Y_cont)
    beta_hsicx_cont <- hsicx_estimator(X_cont, Y_cont, Z_cont, "gaussian")

   beta_ols_nl_cont <- ols_estimator(X_nl_cont, Y_nl_cont)
     beta_iv_nl_cont <- iv_estimator(Z_cont, X_nl_cont, Y_nl_cont)
    beta_hsicx_nl_cont <- hsicx_estimator(X_nl_cont, Y_nl_cont, Z_cont, "gaussian")

    beta_ols_disc <- ols_estimator(X_disc, Y_disc)
    beta_iv_disc <- iv_estimator(Z_disc, X_disc, Y_disc)
    beta_hsicx_disc <- hsicx_estimator(X_disc, Y_disc, Z_disc, "gaussian")

    beta_ols_nl_disc <- ols_estimator(X_nl_disc, Y_nl_disc)
    beta_iv_nl_disc <- iv_estimator(Z_disc, X_nl_disc, Y_nl_disc)
    beta_hsicx_nl_disc <- hsicx_estimator(X_nl_disc, Y_nl_disc, Z_disc, "gaussian")

   
    results[[trial]] <- c(
      beta_ols_cont, beta_iv_cont, beta_hsicx_cont,
         beta_ols_nl_cont, beta_iv_nl_cont, beta_hsicx_nl_cont,
      beta_ols_disc, beta_iv_disc, beta_hsicx_disc,
      beta_ols_nl_disc, beta_iv_nl_disc, beta_hsicx_nl_disc
    )
  }
 

  results_df <- as.data.frame(do.call(rbind, results))
  colnames(results_df) <- c(
    "OLS_Cont", "IV_Cont", "HSICX_Cont",
"OLS_NL_Cont", "IV_NL_Cont", "HSICX_NL_Cont",
    "OLS_Disc", "IV_Disc", "HSICX_Disc",
    "OLS_NL_Disc", "IV_NL_Disc", "HSICX_NL_Disc"
  )
 
  return(results_df)
}

results <- run_simulation(num_trials, n)
print(results)

```
