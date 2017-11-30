library(ensembleEN)
library(MASS)
library(glmnet)
context("Check that the objective function decreases")
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
n <- 100
p <- 200
x_large <- matrix
rho = 0.8
sigma <- (1 - rho) * diag(x = 1, p, p) + rho
x_large <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
true_beta <- c(rep(1, 10), rep(0, p - 10))
y_large <- x_large %*% true_beta + rnorm(n)
x_large_std <- scale(x_large, scale = apply(x_large, 2, function(xj) { sqrt(mean((xj - mean(xj))^2))}))
y_large_cen <- y_large - mean(y_large)
y_large_std <- y_large_cen / sqrt(mean(y_large_cen**2))
# Groups to use for EN in all comparisons
num_groups <- 3
# alpha to use
alphas <- c(0, 1/2, 3/4, 1)
# Number of iterations to perform
num_iter <- c(1, 2, 5, 10, 20)

for (alpha in alphas){
  test_that(paste0("Objective function decreases for p<n, lambda_D=0, ", "alpha=", alpha), {
    fit_glmnet <- glmnet(x_small_std, y_small_std, alpha = alpha)
    lambda_sparsity <- sort(fit_glmnet$lambda)
    length_grid <- length(lambda_sparsity)
    objective_values <- matrix(NA, length(fit_glmnet$lambda), length(num_iter))
    for (i in 1:length(num_iter)){
      fit_phalanx <- Ensemble_EN_Grid(x_small_std, y_small_std, which_lambda = 1, lambdas_grid = lambda_sparsity,
                                                   lambda_fixed = 0, alpha = alpha, num_groups = num_groups,
                                                   tolerance = 0, max_iter = num_iter[i])
      resids_phalanx <- lapply(1:length_grid,
                               function(k, fit_phalanx) { matrix(y_small_std, ncol = num_groups, nrow = length(y_small_std)) - x_small_std %*% fit_phalanx[,,k]}, fit_phalanx)
      
      objective_values[,i] <- sapply(1:length(fit_glmnet$lambda), 
                                       function(k, current_res, beta, lambdas_sparsity, lambdas_diversity, alpha){
                                         Ensemble_EN_Objective(resids_phalanx[[k]], fit_phalanx[,,k], lambda_sparsity[k], 0, alpha)
                                         },
                                         resids_phalanx, fit_phalanx, sort(fit_glmnet$lambda), 0, alpha)
     
    }
    expect_true(all(!apply(objective_values, 2, is.unsorted)))
  })
  
  
  test_that(paste0("Objective function decreases for p>n, lambda_D=0, ", "alpha=", alpha), {
    fit_glmnet <- glmnet(x_large_std, y_large_std, alpha = alpha)
    lambda_sparsity <- sort(fit_glmnet$lambda)
    length_grid <- length(lambda_sparsity)
    objective_values <- matrix(NA, length(fit_glmnet$lambda), length(num_iter))
    for (i in 1:length(num_iter)){
      fit_phalanx <- Ensemble_EN_Grid(x_large_std, y_large_std, which_lambda = 1, lambdas_grid = lambda_sparsity,
                                                   lambda_fixed = 0, alpha = alpha, num_groups = num_groups,
                                                   tolerance = 0, max_iter = num_iter[i])
      resids_phalanx <- lapply(1:length_grid,
                               function(k, fit_phalanx) { matrix(y_large_std, ncol = num_groups, nrow = length(y_large_std)) - x_large_std %*% fit_phalanx[,,k]}, fit_phalanx)
      
      objective_values[,i] <- sapply(1:length(fit_glmnet$lambda), 
                                     function(k, current_res, beta, lambdas_sparsity, lambdas_diversity, alpha){
                                       Ensemble_EN_Objective(resids_phalanx[[k]], fit_phalanx[,,k], lambda_sparsity[k], 0, alpha)
                                     },
                                     resids_phalanx, fit_phalanx, sort(fit_glmnet$lambda), 0, alpha)
      
    }
    expect_true(all(!apply(objective_values, 2, is.unsorted)))
  })
  
  
  test_that(paste0("Objective function decreases for p<n, lambda_D!=0, ", "alpha=", alpha), {
    fit_glmnet <- glmnet(x_large_std, y_large_std, alpha = alpha)
    lambda_sparsity <- sort(fit_glmnet$lambda)
    length_grid <- length(lambda_sparsity)
    objective_values <- matrix(NA, length(fit_glmnet$lambda), length(num_iter))
    for (i in 1:length(num_iter)){
      fit_phalanx <- Ensemble_EN_Grid(x_large_std, y_large_std, which_lambda = 1, lambdas_grid = lambda_sparsity,
                                                   lambda_fixed = 0.1, alpha = alpha, num_groups = num_groups,
                                                   tolerance = 0, max_iter = num_iter[i])
      resids_phalanx <- lapply(1:length_grid,
                               function(k, fit_phalanx) { matrix(y_large_std, ncol = num_groups, nrow = length(y_large_std)) - x_large_std %*% fit_phalanx[,,k]}, fit_phalanx)
      
      objective_values[,i] <- sapply(1:length(fit_glmnet$lambda), 
                                     function(k, current_res, beta, lambdas_sparsity, lambdas_diversity, alpha){
                                       Ensemble_EN_Objective(resids_phalanx[[k]], fit_phalanx[,,k], lambda_sparsity[k], 0, alpha)
                                     },
                                     resids_phalanx, fit_phalanx, sort(fit_glmnet$lambda), 0, alpha)
      
    }
    expect_true(all(!apply(objective_values, 2, is.unsorted)))
  })
  
  
  test_that(paste0("Objective function decreases for p>n, lambda_D!=0, ", "alpha=", alpha), {
    fit_glmnet <- glmnet(x_large_std, y_large_std, alpha = alpha)
    lambda_sparsity <- sort(fit_glmnet$lambda)
    length_grid <- length(lambda_sparsity)
    objective_values <- matrix(NA, length(fit_glmnet$lambda), length(num_iter))
    for (i in 1:length(num_iter)){
      fit_phalanx <- Ensemble_EN_Grid(x_large_std, y_large_std, which_lambda = 1, lambdas_grid = lambda_sparsity,
                                                   lambda_fixed = 0.1, alpha = alpha, num_groups = num_groups,
                                                   tolerance = 0, max_iter = num_iter[i])
      resids_phalanx <- lapply(1:length_grid,
                               function(k, fit_phalanx) { matrix(y_large_std, ncol = num_groups, nrow = length(y_large_std)) - x_large_std %*% fit_phalanx[,,k]}, fit_phalanx)
      
      objective_values[,i] <- sapply(1:length(fit_glmnet$lambda), 
                                     function(k, current_res, beta, lambdas_sparsity, lambdas_diversity, alpha){
                                       Ensemble_EN_Objective(resids_phalanx[[k]], fit_phalanx[,,k], lambda_sparsity[k], 0, alpha)
                                     },
                                     resids_phalanx, fit_phalanx, sort(fit_glmnet$lambda), 0, alpha)
      
    }
    expect_true(all(!apply(objective_values, 2, is.unsorted)))
  })
}