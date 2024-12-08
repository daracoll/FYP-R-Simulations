---
title: "FYP example"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r}
set.seed(348)
true_alpha <- 1
true_beta <- 4
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
Y <- true_alpha * X + true_beta * X^2 + U + epsilon_Y


data <- data.frame(Z = Z, X = X, Y = Y)


X <- data$X 
Y <- data$Y 
Z <- data$Z  


A11 <- sum(Z * X)
A12 <- sum(Z * X^2)
A21 <- sum(Z^2 * X)
A22 <- sum(Z^2 * X^2)

b1 <- sum(Z * Y)
b2 <- sum(Z^2 * Y)


A <- matrix(c(A11, A12, A21, A22), nrow = 2, byrow = TRUE)
b <- c(b1, b2)


params <- solve(A, b)
estimated_alpha <- params[1]
estimated_beta <- params[2]


cat("Estimated alpha:", estimated_alpha, "\n")
cat("Estimated beta:", estimated_beta, "\n")
cat("True alpha:", true_alpha, "\n")
cat("True beta:", true_beta, "\n")

```
