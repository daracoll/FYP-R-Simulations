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
    mse_list[trial, ] <- (params - true_coefficients[1:p])^2  # Standard MSE calculation
  }
  
  if (failed) break
  
  mse_mean <- colMeans(mse_list, na.rm = TRUE)  # Compute standard MSE
  mse_table[[p]] <- mse_mean
  
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

cat("\nTrue Coefficients:\n")
print(true_coefficients)
