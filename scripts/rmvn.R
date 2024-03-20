# Set dimension
dim <- 100

# Create mean vector (1000 zeros)
mean_vector <- rep(0, dim)

# Generate a random matrix
set.seed(123)  # for reproducibility
random_matrix <- matrix(rnorm(dim * dim), nrow = dim)

# Create covariance matrix by multiplying the matrix by its transpose
# This ensures the covariance matrix is symmetric and positive definite
cov_matrix <- random_matrix %*% t(random_matrix)


library(MASS)
library(mvtnorm)
library(microbenchmark)

# Function for generating draws using MASS::mvrnorm
generate_draws_mvrnorm <- function() {
  mvrnorm(n = 1000, mu = mean_vector, Sigma = cov_matrix)
}

# Function for generating draws using mvtnorm::rmvnorm
generate_draws_rmvnorm <- function() {
  rmvnorm(n = 1000, mean = mean_vector, sigma = cov_matrix, method = "chol", checkSymmetry = FALSE)
}

# Function for generating draws using Cholesky Decomposition
generate_draws_chol <- function() {
  chol_decomp <- chol(cov_matrix)
  z <- matrix(rnorm(dim * 1000), 1000, dim)
  t(z %*% chol_decomp)
}

# Benchmarking
benchmark_results <- microbenchmark(
  mvrnorm = generate_draws_mvrnorm(),
  rmvnorm = generate_draws_rmvnorm(),
  chol = generate_draws_chol(),
  times = 10  # Number of times to repeat the test
)

print(benchmark_results)
