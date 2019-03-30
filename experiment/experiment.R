rm(list=ls())
set.seed(310)
dat <- abs(matrix(rexp(40, 1), nrow = 10, ncol = 4))
res <- initialization(dat, max_val = -100, max_iter = 5)

########

k = 2
family = "exponential"
extra_weights = rep(1, nrow(dat))
max_val = -100
max_iter = 10
tol = 1e-3
verbose = F

stopifnot(length(extra_weights) == nrow(dat))
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
pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)
iter <- 1
new_obj <- .evaluate_objective_mat(dat, pred_mat)
old_obj <- Inf

# while(abs(new_obj - old_obj) > tol & iter < max_iter){
#   print(iter)
#   old_obj <- new_obj
#   gradient_mat <- .gradient_mat(dat, pred_mat)
#   new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = k,
#                                      max_val = max_val, direction = direction)
#
#   new_obj <- .evaluate_objective_mat(dat, new_mat)
#   pred_mat <- new_mat
#   iter <- iter + 1
# }
old_obj <- new_obj
gradient_mat <- .gradient_mat(dat, pred_mat)
stepsize_init = 100
stepdown_factor = 2
max_iter = 20

stepsize <- stepsize_init
init_obj <- .evaluate_objective_mat(dat, pred_mat)
iter <- 1

res <- .svd_projection(pred_mat - stepsize*gradient_mat, k = k, factors = T)
# res <- .ensure_feasibility(res$u_mat, res$v_mat, direction = direction,
#                            max_val = max_val, verbose = F)
# new_mat <- res$u_mat %*% t(res$v_mat)
#
# new_obj <- .evaluate_objective_mat(dat, new_mat)

################

u_mat = res$u_mat
v_mat = res$v_mat
range(u_mat %*% t(v_mat))
verbose = T
for(j in 1:nrow(v_mat)){
  res <- .projection_l1(v_mat[j,], u_mat, direction = direction, other_bound = max_val)
  if(attr(res, "status") != 0) break()
  v_mat[j,] <- as.numeric(res)
}
range(u_mat %*% t(v_mat))

if(attr(res, "status") != 0){
  if(verbose) warning("Had to use some janky fixes")
  if(length(which(u_mat < 0)) > length(which(u_mat > 0))) {
    u_mat <- -u_mat; v_mat <- -v_mat
  }
  if(direction == ">=") u_mat <- pmax(u_mat, 1e-3) else u_mat <- pmin(u_mat, -1e-3)

  for(j in 1:nrow(v_mat)){
    v_mat[j,] <- .projection_l1(v_mat[j,], u_mat, direction = direction,
                                other_bound = max_val, tol = 0.1)
  }
}
range(u_mat %*% t(v_mat))

