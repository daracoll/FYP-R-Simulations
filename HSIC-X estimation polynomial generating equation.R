
```{r}
library(kernlab)

# Function to create polynomial features up to degree k (no intercept)
generate_polynomial_features <- function(X, degree) {
  n <- nrow(X)
  poly_X <- NULL  # No intercept term
  
  for (d in 1:degree) {
    poly_X <- cbind(poly_X, X^d)
  }
  
  return(poly_X)
}

# HSIC function with normalisation
compute_HSIC <- function(X, Y, sigma = 1) {
  m <- nrow(X)
  
  # Normalisation
  X <- scale(X)
  Y <- scale(Y)
  
  # Compute RBF kernels
  Kx <- exp(-as.matrix(dist(X))^2 / (2 * sigma^2))
  Ky <- exp(-as.matrix(dist(Y))^2 / (2 * sigma^2))
  
  H <- diag(m) - (1/m) * matrix(1, m, m)  # Centering matrix
  
  HSIC <- sum(H %*% Kx %*% H %*% Ky) / (m^2)
  return(HSIC)
}

# SGD to Minimise HSIC with Polynomial Basis (No Intercept)
sgd_HSIC_polynomial <- function(X, Y, Z, degree = 2, lr = 0.001, epochs = 1000, batch_size = 64) {
  m <- nrow(X)
  
  # Generate polynomial features (no intercept)
  poly_X <- generate_polynomial_features(X, degree)
  p <- ncol(poly_X)  # Number of polynomial coefficients
  
  # Try OLS initialization, fall back to zero vector if error occurs
  theta <- tryCatch({
    solve(t(poly_X) %*% poly_X) %*% (t(poly_X) %*% Y)  # OLS estimation
  }, error = function(e) {
    rep(0, p)  # Return a zero vector if OLS fails
  })
  
  # Adam Optimizer variables
  beta1 <- 0.9
  beta2 <- 0.999
  epsilon <- 1e-8
  m_t <- v_t <- rep(0, p)
  
  for (epoch in 1:epochs) {
    idx <- sample(1:m, batch_size, replace = TRUE)
    X_batch <- poly_X[idx, , drop = FALSE]
    Y_batch <- Y[idx, , drop = FALSE]
    Z_batch <- Z[idx, , drop = FALSE]
    
    # Compute residuals
    residuals <- Y_batch - X_batch %*% theta
    
    # Compute HSIC gradient for each parameter
    grad <- numeric(p)
    eps <- 1e-6  # Stability in finite difference approximation
    
    for (j in 1:p) {
      theta_plus <- theta
      theta_minus <- theta
      theta_plus[j] <- theta[j] + eps
      theta_minus[j] <- theta[j] - eps
      
      res_plus <- Y_batch - X_batch %*% theta_plus
      res_minus <- Y_batch - X_batch %*% theta_minus
      
      HSIC_plus <- compute_HSIC(res_plus, Z_batch)
      HSIC_minus <- compute_HSIC(res_minus, Z_batch)
      
      grad[j] <- (HSIC_plus - HSIC_minus) / (2 * eps)
    }
    
    # Adam update
    m_t <- beta1 * m_t + (1 - beta1) * grad
    v_t <- beta2 * v_t + (1 - beta2) * (grad^2)
    m_hat <- m_t / (1 - beta1^epoch)
    v_hat <- v_t / (1 - beta2^epoch)
    
    theta <- theta - lr * m_hat / (sqrt(v_hat) + epsilon)
    
    if (epoch %% 100 == 0) {
      cat("Epoch:", epoch, "HSIC:", compute_HSIC(residuals, Z_batch), "\n")
    }
  }
  
  return(theta)
}

# Example usage with synthetic polynomial data (No intercept, degree 15, integer coefficients)
set.seed(42)
n <- 200
p <- 1  # Single feature for polynomial basis
degree <- 12  # Degree of polynomial

# Generate integer coefficients between -5 and 5 for the polynomial terms
true_theta <- sample(-5:5, degree, replace = TRUE)

# Generate polynomial data based on true_theta
X <- matrix(rnorm(n * p), ncol = p)
Y <- rep(0, n)

for (i in 1:degree) {
  Y <- Y + true_theta[i] * X^i
}

# Add some noise
Y <- Y + rnorm(n)

# Instrumental variable
Z <- matrix(rnorm(n), ncol = 1)

# Run SGD with HSIC
theta_est <- sgd_HSIC_polynomial(X, Y, Z, degree)
print(theta_est)
true_theta
```
