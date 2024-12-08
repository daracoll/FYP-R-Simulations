```{r}
set.seed(123)

sigma_Z <- 1
sigma_U <- 1
sigma_X <- 1
sigma_Y <- 1
n <- 10000

epsilon_Z <- rnorm(n, mean = 0, sd = sigma_Z)
U <- rnorm(n, mean = 0, sd = sigma_U)
epsilon_X <- rnorm(n, mean = 0, sd = sigma_X)
epsilon_Y <- rnorm(n, mean = 0, sd = sigma_Y)

Z <- epsilon_Z
X <- Z^3 + U + epsilon_X

failed_p <- c()
inv_condition_numbers <- c()
results <- list()

for (p in 1:25) {
  true_coefficients <- runif(p, min = 0.5, max = 5)
  Y <- rowSums(sapply(1:p, function(i) true_coefficients[i] * X^i)) + U + epsilon_Y
  
  A <- matrix(0, nrow = p, ncol = p)
  b <- numeric(p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      A[i, j] <- sum(Z^i * X^j)
    }
    b[i] <- sum(Z^i * Y)
  }
  
  inv_cond_num <- tryCatch({
    cond_num <- kappa(A)
    inv_cond_num <- 1 / cond_num
    params <- solve(A, b)
    results[[as.character(p)]] <- data.frame(
      "True Coefficients" = true_coefficients,
      "Estimated Coefficients" = params
    )
    inv_cond_num
  }, error = function(e) {
    failed_p <<- c(failed_p, p)
    inv_cond_num <- 1 / kappa(A)
    inv_condition_numbers <<- c(inv_condition_numbers, inv_cond_num)
    NULL
  })
}

for (p in names(results)) {
  cat("\nFor p =", p, ":\n")
  print(results[[p]])
}

if (length(failed_p) > 0) {
  plot(
    failed_p, inv_condition_numbers,
    type = "b", pch = 19, col = "red",
    xlab = "Value of p",
    ylab = "Inverse Condition Number",
    main = "Matrix Singularity vs. p"
  )
} else {
  cat("No failures occurred in the range of p.")
}

```

