```{r}
double_factorial <- function(n) {
  if (n <= 0) return(1)
  prod(seq(n, 1, by = -2))
}
sigmaZ <-1
sigmaU <-1
sigmaX<-1
# Input the dimension of the matrix
dim <- 100

A <- matrix(0, nrow = dim, ncol = dim)

for (i in 1:dim) {
  for (j in 1:dim) {
    if ((i + j) %% 2 == 1) {
      # When i + j is odd, A[i, j] = 0 (default value already set)
      next
    }

    if (i %% 2 == 1 && j %% 2 == 1) {
      # Both i and j are odd
      value <- 0
      for (alpha in seq(1, j, by = 2)) {
        for (beta in seq(0, max(0, j - alpha), by = 2)) {  # Ensure valid sequence
          gamma <- j - alpha - beta
          {
            term <- factorial(j) / (factorial(alpha) * factorial(beta) * factorial(gamma)) *
                    double_factorial(i + alpha - 1) *
                    double_factorial(beta - 1) *
                    double_factorial(gamma - 1)*sigmaZ^(i+alpha)*sigmaU^beta*sigmaX^gamma
            value <- value + term
          }
        }
      }
      A[i, j] <- value

    } else if (i %% 2 == 0 && j %% 2 == 0) {
      # Both i and j are even
      value <- 0
      for (alpha in seq(0, j, by = 2)) {
        for (beta in seq(0, max(0, j - alpha), by = 2)) {  # Ensure valid sequence
          gamma <- j - alpha - beta
           {
            term <- factorial(j) / (factorial(alpha) * factorial(beta) * factorial(gamma)) *
                    double_factorial(i + alpha - 1) *
                    double_factorial(beta - 1) *
                    double_factorial(gamma - 1)*sigmaZ^(i+alpha)*sigmaU^beta*sigmaX^gamma
            value <- value + term
          }
        }
      }
      value <- value - double_factorial(i - 1) * double_factorial(j - 1)*sigmaZ^i*sigmaX^j
      A[i, j] <- value
    }
  }
}



A

```


```{r}

compute_inverse_condition_number <- function(matrix) {
  
  return(1 / kappa(matrix))
}


max_p <- dim
inverse_condition_numbers <- numeric(max_p)


for (p in 1:max_p) {
  submatrix <- A[1:p, 1:p]
  inverse_condition_numbers[p] <- compute_inverse_condition_number(submatrix)
}


plot(1:max_p, inverse_condition_numbers, type = "o", col = "blue",
     xlab = "p (Size of Top-Left Submatrix)",
     ylab = "Inverse Condition Number",
     main = "Inverse Condition Number of Top-Left p x p Submatrix vs p")


inverse_condition_numbers

```
