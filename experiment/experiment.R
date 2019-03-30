rm(list=ls())
set.seed(10)
dat <- abs(matrix(rexp(40, 1), nrow = 10, ncol = 4))
# res <- initialization(dat, max_val = -100, max_iter = 5)

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

while(abs(new_obj - old_obj) > tol & iter < max_iter){
  print(iter)
  old_obj <- new_obj
  gradient_mat <- .gradient_mat(dat, pred_mat)
  new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = k,
                                     max_val = max_val, direction = direction)

  new_obj <- .evaluate_objective_mat(dat, new_mat)
  pred_mat <- new_mat
  iter <- iter + 1
}

############

old_obj <- new_obj
gradient_mat <- .gradient_mat(dat, pred_mat)
stepsize_init = 100
stepdown_factor = 2
max_iter = 20

stepsize <- stepsize_init
init_obj <- .evaluate_objective_mat(dat, pred_mat)
iter <- 1

while(TRUE){
  print(iter)
  res <- .svd_projection(pred_mat - stepsize*gradient_mat, k = k, factors = T)
  res <- .ensure_feasibility(res$u_mat, res$v_mat, direction = direction,
                             max_val = max_val, verbose = F)
  new_mat <- res$u_mat %*% t(res$v_mat)

  new_obj <- .evaluate_objective_mat(dat, new_mat)

  print(paste0(init_obj, " : ", new_obj))

  if(new_obj < init_obj) break() else stepsize <- stepsize/stepdown_factor
  iter <- iter + 1
  if(iter > max_iter) stop("Adaptive gradient initialization failed")
}

################

res <- .svd_projection(pred_mat - stepsize*gradient_mat, k = k, factors = T)
u_mat = res$u_mat
v_mat = res$v_mat

j=1
res <- .projection_l1(v_mat[j,], u_mat, direction = direction, other_bound = max_val)

current_vec = v_mat[j,]
other_mat = u_mat
tol = 0.001
other_bound = max_val

stopifnot(ncol(other_mat) == length(current_vec))
stopifnot(is.na(other_bound) || (direction == "<=" & other_bound < 0) ||
            (direction == ">=" & other_bound > 0))

k <- length(current_vec)
d <- nrow(other_mat)
other_direction <- ifelse(direction == "<=", ">=", "<=")

objective_in <- c(rep(1,k), rep(0, k))
constr_mat <- cbind(matrix(0, nrow = nrow(other_mat),ncol = k), other_mat)
constr_mat <- rbind(constr_mat, cbind(diag(k), diag(k)))
constr_mat <- rbind(constr_mat, cbind(diag(k), -diag(k)))

if(direction == "<="){
  constr_ub <- rep(-tol, d)
  if(all(is.na(other_bound))) constr_lb <- rep(-Inf, d) else constr_lb <- rep(other_bound, d)
} else {
  constr_lb <- rep(tol, d)
  if(all(is.na(other_bound))) constr_ub <- rep(Inf, d) else constr_ub <- rep(other_bound, d)
}
constr_lb <- c(constr_lb, current_vec, -current_vec)
constr_ub <- c(constr_ub, rep(Inf, 2*k))

if(all(is.na(other_bound))){
  var_ub <- rep(Inf, 2*k); var_lb <- c(rep(0, k), rep(-Inf, k))
} else {
  var_ub <- c(rep(Inf, k), rep(abs(other_bound), k))
  var_lb <- c(rep(0, k), rep(-abs(other_bound), k))
}

res <- clplite::clp_solve(objective_in, constr_mat, constr_lb, constr_ub, var_lb, var_ub, max = F)

