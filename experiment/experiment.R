rm(-list=ls())
set.seed(10)
true_val <- 1/2
u_mat <- matrix(true_val, nrow = 100, ncol = 1)
v_mat <- -matrix(true_val, nrow = 100, ncol = 1)
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat
class(dat) <- c("gaussian", class(dat)[length(class(dat))])

for(i in 1:nrow(u_mat)){
  for(j in 1:nrow(v_mat)){
    dat[i,j] <- abs(stats::rnorm(1, mean = -1/pred_mat[i,j], sd = -1/(2*pred_mat[i,j])))
  }
}

k = 2
family = "gaussian"
max_val = -1000
max_iter = 10
tol = 1e-3
verbose = F

direction <- .dictate_direction(family)

# initialize
dat <- .matrix_completion(dat, k = k)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

# projected gradient descent
# pred_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
#                                         direction = direction,
#                                         max_iter = max_iter,
#                                         tol = tol)

n <- nrow(dat); d <- ncol(dat)
# pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)

min_val <- min(dat[which(dat > 0)])
dat[which(dat <= 0)] <- min_val/2
pred_mat <- .mean_transformation(dat, family)
direction <- .dictate_direction(family)

class(pred_mat) <- "matrix" #bookeeping purposes
# .nonnegative_matrix_factorization(pred_mat, k = k, direction = direction,
#                                   max_val = max_val)

mat <- pred_mat
if(!is.na(max_val)) stopifnot((direction == "<=" & max_val < 0) | (direction == ">=" & max_val > 0))

# corner cases
if(direction == "<=" & all(mat > 0)) return(NA)
if(direction == ">=" & all(mat < 0)) return(NA)

# prepare the matrix
if(direction == "<=") {
  mat <- -mat
  max_val <- -max_val
}

min_val <- stats::quantile(mat[mat > 0], probs = 0.01)
stopifnot(min_val > 0)
mat[mat < 0] <- min_val
if(!is.na(max_val)) mat[mat > max_val] <- max_val

# perform the nonnegative matrix factorization
stopifnot(all(mat > 0))
res <- NMF::nmf(mat, rank = k)

w_mat <- res@fit@W
h_mat <- res@fit@H
new_mat <- w_mat %*% h_mat

# enforce the max constraint by adjusting one matrix
idx <- unique(which(new_mat > max_val, arr.ind = T)[,2])
if(length(idx) > 0){
  for(j in idx){
    ratio <- max_val/max(new_mat[,j])
    stopifnot(ratio <= 1)
    h_mat[,j] <- h_mat[,j]*ratio
  }

  new_mat <- w_mat %*% h_mat
  stopifnot(all(new_mat <= max_val + 1e-3))
}

stopifnot(new_mat >= 0)
new_mat[new_mat < tol] <- tol
stopifnot(new_mat > 0)

# return the matrix
if(direction == "<=") new_mat <- -new_mat

new_mat
