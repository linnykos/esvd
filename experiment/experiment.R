rm(list=ls())
set.seed(10)
dat <- matrix(abs(rnorm(25, 2, 1)), nrow = 5, ncol = 5)
class(dat) <- c("gaussian", class(dat)[length(class(dat))])
# init <- initialization(dat, family = "gaussian", max_val = 100)

#########

family = "gaussian"
max_val = 100
k = 2
max_iter = 10
tol = 1e-3
verbose = F

direction <- .dictate_direction(family)
if(!is.na(direction) & !is.na(max_val)){
  stopifnot((direction == ">=" & max_val > 0) | (direction == "<=" & max_val < 0))
}

# initialize
dat <- .matrix_completion(dat, k = k)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

# projected gradient descent
pred_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
                                        direction = direction,
                                        max_iter = max_iter,
                                        tol = tol)

#################

n <- nrow(dat); d <- ncol(dat)
# pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)
# iter <- 1
# new_obj <- .evaluate_objective_mat(dat, pred_mat)
# old_obj <- Inf
#
# while(abs(new_obj - old_obj) > tol & iter < max_iter){
#   old_obj <- new_obj
#   gradient_mat <- .gradient_mat(dat, pred_mat)
#   new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = k,
#                                      max_val = max_val, direction = direction)
#
#   new_obj <- .evaluate_objective_mat(dat, new_mat)
#   pred_mat <- new_mat
#   iter <- iter + 1
# }

###########

min_val <- min(dat[which(dat > 0)])
dat[which(dat <= 0)] <- min_val/2
pred_mat <- .mean_transformation(dat, family)
direction <- .dictate_direction(family)

class(pred_mat) <- "matrix" #bookeeping purposes
# .project_rank_feasibility(pred_mat, k = k, direction = direction,
#                           max_val = max_val)$matrix

###############

mat <- pred_mat
stopifnot(!is.na(max_val) | !is.na(direction))
if(!is.na(max_val) & !is.na(direction)) stopifnot((direction == "<=" & max_val < 0) | (direction == ">=" & max_val > 0))

iter <- 1
tol <- ifelse(direction == "<=", -1, 1)

while(iter < max_iter){
  res <- .svd_projection(mat, k = k, factors = T)
  mat <- res$u_mat %*% t(res$v_mat)

  if(is.na(direction)){
    # threshold to be within abs(max_val)
    max_val <- abs(max_val)

    idx <- which(abs(mat) >= max_val)
    val <- mat[idx]
    mat[idx] <- sign(val)*max_val

  } else if (direction == "<=") {
    if(all(mat < 0) && (is.na(max_val) || all(mat > max_val))) return(list(matrix = mat, iter = iter))

    if(any(mat < 0)) tol <- min(tol, stats::quantile(mat[mat < 0], probs = 0.95))
    stopifnot(tol < 0)
    mat[mat > 0] <- tol
    if(!is.na(max_val)) mat[mat < max_val] <- max_val

  } else{
    if(all(mat > 0) && (is.na(max_val) || all(mat < max_val))) return(list(matrix = mat, iter = iter))

    if(any(mat > 0)) tol <- max(tol, stats::quantile(mat[mat > 0], probs = 0.05))
    stopifnot(tol > 0)
    mat[mat < 0] <- tol
    if(!is.na(max_val)) mat[mat > max_val] <- max_val
  }

  iter <- iter + 1
}
