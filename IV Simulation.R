```{r}
set.seed(123)

n <- 10000

beta_xy <- 2
beta_ux <- 3
beta_uy <- 4
beta_zx <- 1.5

sd_u <- 1
sd_x <- 1
sd_y <- 1
sd_z <- 1

U <- rnorm(n, mean = 0, sd = sd_u)
Z <- rnorm(n, mean = 0, sd = sd_z)

X <- beta_zx * Z + beta_ux * U + rnorm(n, mean = 0, sd = sd_x)
Y <- beta_xy * X + beta_uy * U + rnorm(n, mean = 0, sd = sd_y)

naive_model <- lm(Y ~ X)
summary(naive_model)

first_stage <- lm(X ~ Z)
X_hat <- predict(first_stage)

iv_model <- lm(Y ~ X_hat)
summary(iv_model)

cat("\nTrue coefficient for X on Y:", beta_xy, "\n")
cat("Naive estimate (without instrument):", coef(naive_model)["X"], "\n")
cat("IV estimate (with Z as instrument):", coef(iv_model)["X_hat"], "\n")

```
