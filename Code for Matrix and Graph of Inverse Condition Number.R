```{r}
double_factorial <- function(n) {
  if (n <= 0) return(1)
  prod(seq(n, 1, by = -2))
}

sigmaZ <- 1
sigmaU <- 1
sigmaX <- 1
dim <- 3

A <- matrix(0, nrow = dim, ncol = dim)

for (i in 1:dim) {
  for (j in 1:dim) {
    if ((i + j) %% 2 == 1) {
      next
    }

    if (i %% 2 == 1 && j %% 2 == 1) {
      value <- 0
      for (alpha in seq(1, j, by = 2)) {
        for (beta in seq(0, max(0, j - alpha), by = 2)) {
          gamma <- j - alpha - beta
          term <- factorial(j) / (factorial(alpha) * factorial(beta) * factorial(gamma)) *
                  double_factorial(i + alpha - 1) *
                  double_factorial(beta - 1) *
                  double_factorial(gamma - 1) * sigmaZ^(i + alpha) * sigmaU^beta * sigmaX^gamma
          value <- value + term
        }
      }
      A[i, j] <- value

    } else if (i %% 2 == 0 && j %% 2 == 0) {
      value <- 0
      for (alpha in seq(0, j, by = 2)) {
        for (beta in seq(0, max(0, j - alpha), by = 2)) {
          gamma <- j - alpha - beta
          term <- factorial(j) / (factorial(alpha) * factorial(beta) * factorial(gamma)) *
                  double_factorial(beta - 1) *
                  double_factorial(gamma - 1) * sigmaZ^(i + alpha) * sigmaU^beta * sigmaX^gamma * 
                  (double_factorial(i + alpha - 1) - (double_factorial(alpha - 1) * double_factorial(i - 1)))
          value <- value + term
        }
      }
      A[i, j] <- value
    }
  }
}

A

```
