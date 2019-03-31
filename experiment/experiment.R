rm(list=ls())
set.seed(2)
dat <- abs(matrix(rnorm(100), 10, 10))
class(dat) <- c("gaussian", class(dat)[length(class(dat))])
direction <- .dictate_direction(class(dat)[1])

pred_mat <- abs(matrix(rnorm(100), 10, 10))
gradient_mat <-  .gradient_mat(dat, pred_mat)

# new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = 2,
#                                    direction = direction)

k = 2
max_val = NA
stepsize_init = 100
stepdown_factor = 2
max_iter = 20
stepsize <- stepsize_init
init_obj <- .evaluate_objective_mat(dat, pred_mat, ...)
iter <- 1

# while(TRUE){
#   res <- .svd_projection(pred_mat - stepsize*gradient_mat, k = k, factors = T)
#   res <- .ensure_feasibility(res$u_mat, res$v_mat, direction = direction,
#                              max_val = max_val, verbose = F)
#   new_mat <- res$u_mat %*% t(res$v_mat)
#
#   new_obj <- .evaluate_objective_mat(dat, new_mat, ...)
#
#   if(new_obj < init_obj) break() else stepsize <- stepsize/stepdown_factor
#   iter <- iter + 1
#   if(iter > max_iter) stop("Adaptive gradient initialization failed")
# }

res <- .svd_projection(pred_mat - stepsize*gradient_mat, k = k, factors = T)
u_mat = res$u_mat
v_mat = res$v_mat
verbose = T

v_mat2 <- v_mat

for(j in 1:nrow(v_mat)){
  res <- .projection_l1(v_mat[j,], u_mat, direction = direction, other_bound = max_val)
  if(attr(res, "status") != 0) break()
  v_mat2[j,] <- as.numeric(res)
}

if(attr(res, "status") != 0){
  if(verbose) warning("Had to use some janky fixes")
  if(length(which(u_mat < 0)) > length(which(u_mat > 0))) {
    u_mat <- -u_mat; v_mat <- -v_mat
  }
  if(direction == ">=") {
    u_mat <- pmin(pmax(u_mat, 1e-3), max_val)
  } else {
    u_mat <- pmax(pmin(u_mat, -1e-3), max_val)
  }

  for(j in 1:nrow(v_mat)){
    v_mat[j,] <- .projection_l1(v_mat[j,], u_mat, direction = direction,
                                other_bound = max_val)
  }
} else {
  v_mat <- v_mat2
}

