context("Test SBM")

## .sbm_projection is correct

test_that(".sbm_projection works (symmetric case)", {
  set.seed(10)
  n <- 100
  b_mat <- matrix(c(0.9, 0.2, 0.2, 0.9), 2, 2)
  cluster_vec <- rep(1:2, each = n/2)

  p_mat <- matrix(NA, n, n)
  mat <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      p_mat[i,j] <- b_mat[cluster_vec[i], cluster_vec[j]]
      p_mat[j,i] <- p_mat[i,j]

      mat[i,j] <- stats::rbinom(1, 1, p = p_mat[i,j])
      mat[j,i] <- mat[i,j]
    }
  }

  res <- .sbm_projection(mat, k = 2)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(mat)))
})

test_that(".sbm_projection works (asymmetric case)", {
  set.seed(10)
  n <- 500; d <- 250
  b_mat <- matrix(c(0.9, 0.2, 0.1, 0.8), 2, 2)
  cluster_row <- rep(1:2, each = n/2)
  cluster_col <- rep(1:2, each = n/2)

  p_mat <- matrix(NA, n, d)
  mat <- matrix(NA, n, d)
  for(i in 1:n){
    for(j in 1:d){
      p_mat[i,j] <- b_mat[cluster_row[i], cluster_col[j]]

      mat[i,j] <- stats::rbinom(1, 1, p = p_mat[i,j])
    }
  }

  res <- .sbm_projection(mat, k = 2)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(mat)))
})

test_that(".sbm_projection keeps an estimate that is within the range of mat", {
  load("../assets/sbm1.RData")
  res <- .sbm_projection(mat, k = 2)

  expect_true(min(res) > 0)
  expect_true(max(res) <= max(mat))
})


test_that(".sbm_projection optimizes its respective model", {
  trials <- 20
  n <- 100; d <- 80; k <- 2

  dat_list <- lapply(1:trials, function(i){
    set.seed(i)
    b_mat <- matrix(sample(10:100, 4), 2, 2)

    cluster_row <- sample(1:2, size = n, replace = T)
    cluster_col <- sample(1:2, size = d, replace = T)

    p_mat <- matrix(NA, n, d)
    mat <- matrix(NA, n, d)
    for(i in 1:n){
      for(j in 1:d){
        p_mat[i,j] <- b_mat[cluster_row[i], cluster_col[j]]

        mat[i,j] <- abs(stats::rnorm(1, p_mat[i,j]))
      }
    }

    list(p_mat = p_mat, mat = mat)
  })

  fit_list <- lapply(1:length(dat_list), function(i){
    set.seed(i)
    .sbm_projection(dat_list[[i]]$mat, k = 2)
  })

  bool_vec <- sapply(1:length(dat_list), function(i){
    min(fit_list[[i]] > 0) & max(fit_list[[i]]) < max(dat_list[[i]]$mat)
  })

  expect_true(all(bool_vec))

  error_mat <- sapply(1:length(dat_list), function(i){
    vec <- sapply(1:length(fit_list), function(j){
      .l2norm(fit_list[[j]] - dat_list[[i]]$p_mat)
    })

    vec/min(vec)
  })

  expect_true(all(diag(error_mat) <= 1+1e-6))
})
