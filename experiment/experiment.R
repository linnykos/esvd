rm(list=ls())
i = 10
set.seed(i)
n <- 20
u_mat <- matrix(abs(rnorm(n)), ncol = 1)
v_mat <- matrix(abs(rnorm(n)), ncol = 1)
pred_mat <- u_mat %*% t(v_mat)

dat <- pred_mat

for(i in 1:n){
  for(j in 1:n){
    dat[i,j] <- abs(stats::rnorm(1, mean = 1/pred_mat[i,j], sd = 1/(2*pred_mat[i,j])))
  }
}

class(dat) <- c("curved_gaussian", class(dat)[length(class(dat))])

# init <- initialization(dat, k = 1, family = "curved_gaussian", max_val = 100)

############

k = 1
family= "curved_gaussian"
max_val = 100
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
# pred_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
#                                         direction = direction,
#                                         max_iter = max_iter,
#                                         tol = tol)

##################

n <- nrow(dat); d <- ncol(dat)
# pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)
# iter <- 1
# new_obj <- .evaluate_objective_mat(dat, pred_mat)
# old_obj <- Inf
#
# while(abs(new_obj - old_obj) > tol & iter < max_iter){
#   old_obj <- new_obj
#   gradient_mat <- .gradient_mat(dat, pred_mat, ...)
#   new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = k,
#                                      max_val = max_val, direction = direction)
#
#   new_obj <- .evaluate_objective_mat(dat, new_mat)
#   pred_mat <- new_mat
#   iter <- iter + 1
# }
#

##########################

min_val <- max(min(dat[which(dat > 0)]), 1e-3)
dat[which(dat <= min_val)] <- min_val/2
pred_mat <- .mean_transformation(dat, family)
direction <- .dictate_direction(family)

plot(sort(as.numeric(pred_ma)))

# .project_rank_feasibility(pred_mat, k = k, direction = direction,
#                           max_val = max_val)$matrix

###############

mat <- pred_mat
stopifnot(!is.na(max_val) | !is.na(direction))
if(!is.na(max_val) & !is.na(direction)) stopifnot((direction == "<=" & max_val < 0) | (direction == ">=" & max_val > 0))

mat_org <- mat
iter <- 1

while(iter < 3){
  res <- .svd_projection(mat, k = k, factors = T)
  mat <- res$u_mat %*% t(res$v_mat)

  if(is.na(direction)){
    if(all(abs(mat) <= max_val+tol)) return(list(matrix = mat, iter = iter))
  } else if (direction == "<=") {
    if(all(mat < 0) && (is.na(max_val) || all(mat > max_val-tol))) return(list(matrix = mat, iter = iter))
  } else {
    if(all(mat > 0) && (is.na(max_val) || all(mat < max_val+tol))) return(list(matrix = mat, iter = iter))
  }

  mat <- .absolute_threshold(mat, direction, max_val)

  iter <- iter + 1
}

# res <- .svd_projection(mat, k = k, factors = T)

#############

res <- RSpectra::svds(mat, k = k)
