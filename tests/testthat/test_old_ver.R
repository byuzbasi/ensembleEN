library(ensembleEN)
library(MASS)
context("Compare main algorithm with saved old output")
# Generate data sets, one with p<n, the other with p>n
set.seed(1)
n <- 100
p <- 20
x_small <- matrix
rho = 0.8
sigma <- (1 - rho) * diag(x = 1, p, p) + rho
x_small <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
true_beta <- c(rep(1, 10), rep(0, p - 10))
y_small <- x_small %*% true_beta + rnorm(n)
x_small_std <- scale(x_small, scale = apply(x_small, 2, function(xj) { sqrt(mean((xj - mean(xj))^2))}))
y_small_cen <- y_small - mean(y_small)
y_small_std <- y_small_cen / sqrt(mean(y_small_cen**2))

set.seed(1)
n <- 50
p <- 100
x_large <- matrix
rho = 0.95
sigma <- (1 - rho) * diag(x = 1, p, p) + rho
x_large <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
true_beta <- c(rep(1, 10), rep(0, p - 10))
y_large <- x_large %*% true_beta + rnorm(n)
x_large_std <- scale(x_large, scale = apply(x_large, 2, function(xj) { sqrt(mean((xj - mean(xj))^2))}))
y_large_cen <- y_large - mean(y_large)
y_large_std <- y_large_cen / sqrt(mean(y_large_cen**2))
# alpha to use
alphas <- c(0.01, 1/2, 3/4, 1)

for(alpha in alphas){
  test_that(paste0("Equality for p<n, alpha = ", alpha), {
    # Load saved fit
    file_name <- paste0('fits_test_small_', alpha, '.Rdata')
    load(file_name)
    # Recompute
    set.seed(1)
    fit_phalanx_new <- ensembleEN(x_small_std, y_small_std, alpha=alpha)
    # Compare all entries of the new and old fit
    diffs <- sapply(names(fit_phalanx), function(k, fit_phalanx, fit_phalanx_new){
      sum((fit_phalanx[[k]] - fit_phalanx_new[[k]])**2)
    }, fit_phalanx, fit_phalanx_new)
    print(diffs)
    expect_true(max(diffs) < 1e-20)
  })
  
  test_that(paste0("Equality for p>n, alpha = ", alpha), {
    # Load saved fit
    file_name <- paste0('fits_test_large_', alpha, '.Rdata')
    load(file_name)
    # Recompute
    set.seed(1)
    fit_phalanx_new <- ensembleEN(x_large_std, y_large_std, alpha=alpha)
    # Compare all entries of the new and old fit
    diffs <- sapply(names(fit_phalanx), function(k, fit_phalanx, fit_phalanx_new){
      sum((fit_phalanx[[k]] - fit_phalanx_new[[k]])**2)
    }, fit_phalanx, fit_phalanx_new)
    expect_true(max(diffs) < 1e-20)
    print(diffs)
  })
}