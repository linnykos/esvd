rm(list=ls())
x <- 23
set.seed(x*10)
dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))

res <- initialization(dat, max_val = -100)
range(res$u_mat %*% t(res$v_mat))

######

family = "exponential"
max_iter = 10
tol = 1e-3
verbose = F
max_val = -100
direction <- .dictate_direction(family)
k = 2

# initialize
dat <- .matrix_completion(dat, k = k)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

# projected gradient descent
pred_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
                                        direction = direction,
                                        max_iter = max_iter,
                                        tol = tol)

##########

n <- nrow(dat); d <- ncol(dat)
pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)

###########

min_val <- min(dat[which(dat > 0)])
dat[which(dat <= 0)] <- min_val/2
pred_mat <- .mean_transformation(dat, family)
direction <- .dictate_direction(family)

class(pred_mat) <- "matrix" #bookeeping purposes
zz <- .nonnegative_matrix_factorization(pred_mat, k = k, direction = direction,
                                  max_val = max_val)

##########

mat <- pred_mat
tol = 1e-3

# corner cases
if(direction == "<=" & all(mat > 0)) return(NA)
if(direction == ">=" & all(mat < 0)) return(NA)

# prepare the matrix
if(direction == "<=") {
  mat <- -mat
  max_val <- -max_val
}

min_val <- quantile(mat[mat > 0], probs = 0.01)
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

