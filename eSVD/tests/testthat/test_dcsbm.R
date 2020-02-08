context("Testing DCSBM")

## .dcsbm_projection is correct

test_that(".dcsbm_projection works (symmetric case)", {
  set.seed(10)
  n <- 100
  theta_vec <- runif(n, min = 0.5, max = 1)
  b_mat <- matrix(c(0.9, 0.2, 0.2, 0.9), 2, 2)
  cluster_vec <- rep(1:2, each = n/2)

  p_mat <- matrix(NA, n, n)
  mat <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      p_mat[i,j] <- theta_vec[i]*theta_vec[j]*b_mat[cluster_vec[i], cluster_vec[j]]
      p_mat[j,i] <- p_mat[i,j]

      mat[i,j] <- stats::rbinom(1, 1, p = p_mat[i,j])
      mat[j,i] <- mat[i,j]
    }
  }

  res <- .dcsbm_projection(mat, k = 2)

  expect_true(is.list(res))
  expect_true(length(res) == 4)
  expect_true(all(dim(res$base_mat) == c(n,n)))
})

test_that(".dcsbm_projection works (asymmetric case)", {
  set.seed(10)
  n <- 500; d <- 250
  theta_row <- runif(n, min = 0.5, max = 1)
  theta_col <- runif(d, min = 0.5, max = 1)
  b_mat <- matrix(c(0.9, 0.2, 0.1, 0.8), 2, 2)
  cluster_row <- rep(1:2, each = n/2)
  cluster_col <- rep(1:2, each = n/2)

  p_mat <- matrix(NA, n, d)
  mat <- matrix(NA, n, d)
  for(i in 1:n){
    for(j in 1:d){
      p_mat[i,j] <- theta_row[i]*theta_col[j]*b_mat[cluster_row[i], cluster_col[j]]

      mat[i,j] <- stats::rbinom(1, 1, p = p_mat[i,j])
    }
  }

  res <- .dcsbm_projection(mat, k = 2)

  expect_true(is.list(res))
  expect_true(length(res) == 4)
  expect_true(all(dim(res$base_mat) == c(n,n)))
})
