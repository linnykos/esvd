rm(list=ls())

loss_gaussian_pop <- function(nat_mat, u_mat, v_mat, scalar = 2, ...){
  n <- nrow(nat_mat); p <- ncol(nat_mat)
  pred_mat <- u_mat %*% t(v_mat)

  val <- sum(-log(-pred_mat) +
        pred_mat*(-1/nat_mat)*scalar^2 +
        pred_mat^2*((1/scalar^2 + 1)*1/nat_mat^2)*scalar^2/2)

  val/p
}

loss_gaussian_samp <- function(dat, u_mat, v_mat, scalar = 2, ...){
  n <- nrow(nat_mat); p <- ncol(nat_mat)
  pred_mat <- u_mat %*% t(v_mat)

  val <- sum(-log(-pred_mat) +
               pred_mat*dat*scalar^2 +
               pred_mat^2*dat^2*scalar^2/2)

  val/p
}

# here, u_mat is always the left
gradient_gaussian_pop <- function(nat_mat, u_mat, v_mat, scalar = 2, ...){
  n <- nrow(u_mat); p <- nrow(v_mat); k <- ncol(u_mat)
  pred_mat <- u_mat %*% t(v_mat)

  grad_mat <- matrix(0, n, k)

  for(i in 1:n){
    grad_mat[i,] <- rowSums(sapply(1:p, function(j){
      (-1/(pred_mat[i,j]) + scalar^2*(-1/nat_mat[i,j]) +
         scalar^2*((1/scalar^2 + 1)*1/nat_mat[i,j]^2)*pred_mat[i,j]) * v_mat[j,]
    }))
  }

  grad_mat/p
}

gradient_gaussian_samp <- function(dat, u_mat, v_mat, scalar = 2, ...){
  n <- nrow(u_mat); p <- nrow(v_mat); k <- ncol(u_mat)
  pred_mat <- u_mat %*% t(v_mat)

  grad_mat <- matrix(0, n, k)

  for(i in 1:n){
    grad_mat[i,] <- rowSums(sapply(1:p, function(j){
      (-1/(pred_mat[i,j]) + scalar^2*dat[i,j] +
         scalar^2*dat[i,j]^2*pred_mat[i,j]) * v_mat[j,]
    }))
  }

  grad_mat/p
}

# here, u_mat is always the left
hessian_gaussian_pop <- function(nat_mat, u_mat, v_mat,  n, p, scalar = 2, ...){
  n <- nrow(u_mat); p <- nrow(v_mat); k <- ncol(u_mat)
  pred_mat <- u_mat %*% t(v_mat)

  hes_mat <- matrix(0, n*k, n*k)

  for(i in 1:n){
    tmp <- matrix(0, k, k)

    for(j in 1:p){
      tmp <- tmp + (1/(pred_mat[i,j]^2) + scalar^2*((1/scalar^2 + 1)*1/nat_mat[i,j]^2)) *
        v_mat[j,] %*% t(v_mat[j,])
    }

    hes_mat[((i-1)*k+1):(i*k), ((i-1)*k+1):(i*k)] <- tmp
  }

  hes_mat/p
}

hessian_gaussian_samp <- function(dat, u_mat, v_mat,  n, p, scalar = 2, ...){
  n <- nrow(u_mat); p <- nrow(v_mat); k <- ncol(u_mat)
  pred_mat <- u_mat %*% t(v_mat)

  hes_mat <- matrix(0, n*k, n*k)

  for(i in 1:n){
    tmp <- matrix(0, k, k)

    for(j in 1:p){
      tmp <- tmp + (1/(pred_mat[i,j]^2) + scalar^2*dat[i,j]^2) *
        v_mat[j,] %*% t(v_mat[j,])
    }

    hes_mat[((i-1)*k+1):(i*k), ((i-1)*k+1):(i*k)] <- tmp
  }

  hes_mat/p
}

generate_pairs <- function(n, p, k){
  # left_mat <- abs(matrix(runif(n*k), n, k))
  # right_mat <- -abs(matrix(runif(p*k), p, k))
  val <- 1 # runif(1)
  left_mat <- abs(matrix(val, n, k))
  right_mat <- -abs(matrix(val, p, k))

  nat_mat <- left_mat %*% t(right_mat)
  svd_res <- svd(nat_mat)
  u_mat <- svd_res$u[,1:k] %*% diag(svd_res$d[1:k])/sqrt(p)
  v_mat <- svd_res$v[,1:k]*sqrt(p)

  dat <- matrix(0, n, p)
  for(i in 1:n){
    for(j in 1:p){
      theta <- -1/nat_mat[i,j]
      dat[i,j] <- rnorm(1, mean = theta, sd = theta/scalar)
    }
  }

  list(u_mat = u_mat, v_mat = v_mat, nat_mat = nat_mat, dat = dat)
}

###########################

set.seed(10)
n <- 20; p <- 500
#n <- 100; p <- 50
#n <- 500; p <- 200
k <- 3; scalar <- 2
res <- generate_pairs(n, p, k)
res2 <- generate_pairs(n, p, k)

nat_mat <- res$nat_mat
dat <- res$dat
u_mat <- res2$u_mat
v_mat <- res2$v_mat

# propose a different set

l_pop <- loss_gaussian_pop(nat_mat, u_mat, v_mat)
g_pop <- gradient_gaussian_pop(nat_mat, u_mat, v_mat)
h_pop <- hessian_gaussian_pop(nat_mat, u_mat, v_mat)

range(eigen(h_pop)$values)

g_samp <- gradient_gaussian_samp(dat, u_mat, v_mat)
.l2norm(g_pop - g_samp)

l_samp <- loss_gaussian_samp(dat, u_mat, v_mat)
abs(l_pop - l_samp)

################################

# check gradient is correct
trials <- 1000
k <- 3; scalar <- 2
n <- 100; p <- 50
set.seed(0)
res <- generate_pairs(n, p, k)
nat_mat <- res$nat_mat

bool <- sapply(1:trials, function(x){
  if(x %% floor(trials/10) == 0) cat('*')
  set.seed(x)
  left1 <- abs(matrix(runif(n*k), n, k))
  left2 <- abs(matrix(runif(n*k), n, k))
  right1 <- -abs(matrix(runif(p*k), p, k))

  grad1 <- gradient_gaussian_pop(nat_mat, left1, right1)
  grad2 <- gradient_gaussian_pop(nat_mat, left2, right1)

  diff1 <- grad1 - grad2
  diff2 <- left1 - left2

  sol <- psych::tr(t(diff1) %*% diff2)

  sol >= 0
})

# ######################################
#
# # check the hessian is correct
# trials <- 1000
# k <- 3; scalar <- 2
# n <- 100; p <- 50
# set.seed(0)
# res <- generate_pairs(n, p, k)
# nat_mat <- res$nat_mat
#
# bool <- sapply(1:trials, function(x){
#   if(x %% floor(trials/10) == 0) cat('*')
#   set.seed(x)
#   left1 <- abs(matrix(runif(n*k), n, k))
#   left2 <- abs(matrix(runif(n*k), n, k))
#   right1 <- -abs(matrix(runif(p*k), p, k))
#
#   grad1 <- gradient_gaussian_pop(nat_mat, left1, right1)
#   grad2 <- gradient_gaussian_pop(nat_mat, left2, right1)
#
#   hess1 <- hessian_gaussian_pop(nat_mat, left1, right1)
#   hess2 <- hessian_gaussian_pop(nat_mat, left2, right1)
#
#   diff1 <- hess1 - hess2
#   diff2 <- grad1 - grad2
#
#   sol <- psych::tr(t(diff1) %*% diff2)
#
#   sol >= 0
# })
