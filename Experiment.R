

```{r}
library(AER)      
library(dHSIC)   
library(dplyr)    
library(lmtest)   
library(sandwich) 
library(kernlab)  

data("CigarettesSW", package = "AER")  
CigarettesSW <- as.data.frame(na.omit(CigarettesSW))  

Y_cig <- CigarettesSW$packs  
X_cig <- CigarettesSW$price  
Z_cig <- CigarettesSW$tax    

first_stage_cig <- lm(X_cig ~ Z_cig, data = CigarettesSW)
summary(first_stage_cig)

first_stage_fstat_cig <- summary(first_stage_cig)$fstatistic[1]
cat("First-Stage F-Statistic (CIGARETTES):", first_stage_fstat_cig, "\n")

iv_model_cig <- ivreg(Y_cig ~ X_cig | Z_cig, data = CigarettesSW)
summary(iv_model_cig)
theta_2sls_cig <- coef(iv_model_cig)  

hsic_grad <- function(x, y, z, theta){
  theta1 <- theta[2] - 0.001
  theta2 <- theta[2] + 0.001
  resid1 <- y - (theta[1] + theta1 * x)
  resid2 <- y - (theta[1] + theta2 * x)
  hsic1 <- dhsic(matrix(resid1, ncol = 1), matrix(z, ncol = 1))$dHSIC
  hsic2 <- dhsic(matrix(resid2, ncol = 1), matrix(z, ncol = 1))$dHSIC
  return((hsic2 - hsic1) / (theta2 - theta1))
}

hsicx <- function(X, Y, Z, t = 100, k = 100, m = 32, gamma = 0.1, alpha = .05) {
  theta_best <- c(theta_2sls_cig[1], theta_2sls_cig[2])  
  pval_best <- 0
  l <- 0  
  pval <- 0  
  
  while (l < t && pval < alpha) {
    if (l == 0) {
      theta <- c(theta_2sls_cig[1], theta_2sls_cig[2])  
    } else {
      theta <- c(rnorm(1, mean = theta_2sls_cig[1], sd = 10), 
                 rnorm(1, mean = theta_2sls_cig[2], sd = 10))
    }
    
    update <- 1  
    while (update > 0.001) {  
      theta_old <- theta[2]
      
      for (i in 1:k) {
        sample_data <- slice_sample(data.frame(X, Y, Z), n = min(m, nrow(CigarettesSW)))  
        theta[2] <- theta[2] - gamma * hsic_grad(sample_data$X, sample_data$Y, sample_data$Z, theta)
      }
      
      update <- sqrt((theta_old - theta[2])^2)  
    }
    
    resid <- matrix(Y - (theta[1] + theta[2] * X), ncol = 1)  
    test <- dhsic.test(resid, matrix(Z, ncol = 1))  
    pval <- test$p.value  
    
    cat("Iteration:", l + 1, "| Intercept:", theta[1], "| Coefficient:", theta[2], "| P-value:", pval, "\n")
    
    if (pval > alpha) {
      cat("Stopping early: P-value exceeded alpha at iteration", l + 1, "\n")
      return(theta)  
    }
    
    if (pval > pval_best) {  
      pval_best <- pval
      theta_best <- theta
    }
    
    l <- l + 1  
  }
  
  return(theta_best)
}

theta_hsicx_cig <- hsicx(X_cig, Y_cig, Z_cig)

cat("2SLS Estimate (CIGARETTES): Intercept:", theta_2sls_cig[1], "Coefficient:", theta_2sls_cig[2], "\n")
cat("HSIC-X Estimate (CIGARETTES): Intercept:", theta_hsicx_cig[1], "Coefficient:", theta_hsicx_cig[2], "\n")

set.seed(123)  
n_boot <- 100
boot_2sls <- numeric(n_boot)

for (i in 1:n_boot) {
  sample_indices <- sample(1:nrow(CigarettesSW), replace = TRUE)
  sample_data <- CigarettesSW[sample_indices, ]
  
  boot_model <- ivreg(packs ~ price | tax, data = sample_data)
  boot_2sls[i] <- coef(boot_model)[2]  
}

ci_2sls <- quantile(boot_2sls, probs = c(0.025, 0.975))

set.seed(123)
boot_hsicx <- numeric(n_boot)

for (i in 1:n_boot) {
  sample_indices <- sample(1:nrow(CigarettesSW), replace = TRUE)
  sample_data <- CigarettesSW[sample_indices, ]
  
  boot_hsicx[i] <- hsicx(sample_data$price, sample_data$packs, sample_data$tax)[2]  
}

ci_hsicx <- quantile(boot_hsicx, probs = c(0.025, 0.975))

cat("2SLS Estimate:", theta_2sls_cig[2], "\n")
cat("2SLS 95% CI:", ci_2sls, "\n")
cat("HSIC-X Estimate:", theta_hsicx_cig[2], "\n")
cat("HSIC-X 95% CI:", ci_hsicx, "\n")


```

