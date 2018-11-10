.fit_factorization <- function(dat, u_mat, v_mat, max_val = NA,
                               family = "exponential",
                               reparameterize = F,
                               extra_weights = rep(1, nrow(dat)),
                               tol = 1e-3, max_iter = 100,
                               verbose = F,
                               cores = NA){
  if(!is.na(cores)) doMC::registerDoMC(cores = cores)
  stopifnot(length(which(dat > 0)) > 0)
  stopifnot(family == "poisson" | all(extra_weights == 1))
  stopifnot(all(extra_weights > 0))
  stopifnot(length(which(dat < 0)) == 0)
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))
  k <- ncol(u_mat)
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

  idx <- which(dat == 0)
  min_val <- min(dat[which(dat > 0)])
  dat[which(dat == 0)] <- min_val/2

  current_obj <- Inf
  next_obj <- .evaluate_objective(dat, u_mat, v_mat, extra_weights = extra_weights)
  obj_vec <- c(next_obj)
  if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))

  while((is.na(tol) | abs(current_obj - next_obj) > tol) & length(obj_vec) < max_iter){
    current_obj <- next_obj

    u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, extra_weights = extra_weights,
                           !is.na(cores))
    v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, extra_weights = extra_weights,
                           !is.na(cores))

    next_obj <- .evaluate_objective(dat, u_mat, v_mat, extra_weights = extra_weights)

    if(reparameterize){
      pred_mat <- u_mat %*% t(v_mat)
      svd_res <- svd(pred_mat)
      if(k == 1) diag_vec <- as.matrix(sqrt(svd_res$d[1:k])) else diag_vec <- sqrt(svd_res$d[1:k])
      u_mat <- svd_res$u[,1:k] %*% diag(diag_vec)
      v_mat <- svd_res$v[,1:k] %*% diag(diag_vec)
    }

    if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))

    obj_vec <- c(obj_vec, next_obj)
  }

  pred_mat <- u_mat %*% t(v_mat)
  svd_res <- svd(pred_mat)

  if(k == 1) diag_vec <- as.matrix(sqrt(svd_res$d[1:k])) else diag_vec <- sqrt(svd_res$d[1:k])
  list(u_mat = svd_res$u[,1:k,drop = F] %*% diag(diag_vec),
       v_mat = svd_res$v[,1:k,drop = F] %*% diag(diag_vec),
       obj_vec = obj_vec)
}

#########

.evaluate_objective <- function (dat, u_mat, v_mat, ...) {
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))

  UseMethod(".evaluate_objective")
}

.evaluate_objective.default <- function(dat, u_mat, v_mat, ...){
  .evaluate_objective.exponential(dat, u_mat, v_mat, ...)
}

.evaluate_objective_single <- function (dat_vec, current_vec, other_mat, ...) {
  stopifnot(!is.matrix(dat_vec))
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  UseMethod(".evaluate_objective_single")
}

.evaluate_objective_single.default <- function(dat_vec, current_vec, other_mat, ...){
  .evaluate_objective_single.exponential(dat_vec, current_vec, other_mat, ...)
}

.gradient_vec <- function (dat_vec, current_vec, other_mat, ...) {
  stopifnot(!is.matrix(dat_vec))
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  UseMethod(".gradient_vec")
}

.gradient_vec.default <- function(dat_vec, current_vec, other_mat, ...){
  .gradient_vec.exponential(dat_vec, current_vec, other_mat, ...)
}

#########

.optimize_mat <- function(dat, current_mat, other_mat, left = T, max_val = NA,
                          extra_weights = rep(1, nrow(dat)), parallelized = F){
  stopifnot(length(class(dat)) == 2)

  stopifnot(ncol(current_mat) == ncol(other_mat))
  if(left) {
    stopifnot(ncol(dat) == nrow(other_mat))
  } else {
    stopifnot(nrow(dat) == nrow(other_mat))
  }

  if(parallelized){
    func <- function(i){
      if(left) {
        dat_vec <- dat[i,]; extra_vec <- rep(extra_weights[i], nrow(other_mat))
      } else {
        dat_vec <- dat[,i]; extra_vec <- extra_weights
      }
      class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
      .optimize_row(dat_vec, current_mat[i,], other_mat, max_val = max_val, extra_weights = extra_vec)
    }

    lis <- foreach::"%dopar%"(foreach::foreach(i = 1:nrow(current_mat)), func(i))
    current_mat <- do.call(rbind, lis)

  } else {
    for(i in 1:nrow(current_mat)){
      if(left) {
        dat_vec <- dat[i,]; extra_vec <- rep(extra_weights[i], nrow(other_mat))
      } else {
        dat_vec <- dat[,i]; extra_vec <- extra_weights
      }
      class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
      if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                                other_mat, max_val = max_val,
                                                                extra_weights = extra_vec)
    }
  }

  current_mat
}

.optimize_row <- function(dat_vec, current_vec, other_mat, max_iter = 100,
                          max_val = NA, extra_weights = rep(1, nrow(other_mat))){
  stopifnot(length(which(!is.na(dat_vec))) > 0)
  stopifnot(length(extra_weights) == nrow(other_mat))

  direction <- .dictate_direction(class(dat_vec)[1])
  current_obj <- Inf
  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat, extra_weights = extra_weights)
  iter <- 1

  while(abs(current_obj - next_obj) > 1e-6 & iter < max_iter){
    current_obj <- next_obj

    grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat, extra_weights = extra_weights)
    step_vec <- .frank_wolfe(grad_vec, other_mat, which(!is.na(dat_vec)),
                             direction = direction, other_bound = max_val)
    step_size <- .binary_search(dat_vec, current_vec, step_vec, other_mat, extra_weights = extra_weights)
    current_vec <- (1-step_size)*current_vec + step_size*step_vec

    next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat, extra_weights = extra_weights)

    iter <- iter + 1
  }

  current_vec
}

.binary_search <- function(dat_vec, current_vec, step_vec, other_mat,
                           max_iter = 100, extra_weights = rep(1, nrow(other_mat))){
  form_current <- function(s){
    stopifnot(0<=s, s<=1)
    (1-s)*current_vec + s*step_vec
  }

  upper <- 1; lower <- 0; mid <- 0.5
  iter <- 1

  upper_val <- .evaluate_objective_single(dat_vec, form_current(upper), other_mat, extra_weights = extra_weights)
  lower_val <- .evaluate_objective_single(dat_vec, form_current(lower), other_mat, extra_weights = extra_weights)

  upper_val_org <- upper_val; lower_val_org <- lower_val
  mid_val <- 2*max(upper_val, lower_val);

  # stage 1: do at least a few iterations first
  while(iter <= 6){
    mid_val <- .evaluate_objective_single(dat_vec, form_current(mid), other_mat, extra_weights = extra_weights)

    if(lower_val <= upper_val){
      upper <- mid; upper_val <- mid_val
    } else {
      lower <- mid; lower_val <- mid_val
    }

    mid <- (upper + lower)/2
    iter <- iter + 1
  }

  # stage 2: ensure that the stepsize actually decreases the obj val
  while(iter <= max_iter & mid_val > min(upper_val_org, lower_val_org)){
    mid_val <- .evaluate_objective_single(dat_vec, form_current(mid), other_mat, extra_weights = extra_weights)

    if(lower_val <= upper_val){
      upper <- mid; upper_val <- mid_val
    } else {
      lower <- mid; lower_val <- mid_val
    }

    mid <- (upper + lower)/2
    iter <- iter + 1
  }

  if(upper_val_org <= mid_val) return(1)
  if(lower_val_org <= mid_val) return(0)

  mid
}

.l2norm <- function(x){sqrt(sum(x^2))}

#######################

#' Frank wolfe linear program
#'
#' Solves the linear program:
#' Data: g_i (vector), i from 1 to k. V (matrix), with d rows and k columns.
#' Variables: y_+1 to y_+k, y_-1 to y_-k
#' Optimization problem: min g_1*(y_+1-y_-1) + ... + g_k(y_+k-y_-k)
#' such that: V %*% (y_+(1:k) - y_-(1:k)) <= tol elementwise (for entire vector of length d)
#'            V %*% (y_+(1:k) - y_-(1:k)) >= other_bound
#'            y_+i, y_-i >= 0 for i from 1 to k
#'
#' @param grad_vec vector
#' @param other_mat matrix
#' @param idx row indices for other_mat
#' @param tol numierc
#' @param direction character
#' @param other_bound numeric
#'
#' @return vector
.frank_wolfe <- function(grad_vec, other_mat, idx = 1:nrow(other_mat),
                         tol = 0.0001, direction = "<=", other_bound = NA){

  k <- length(grad_vec)
  other_mat <- other_mat[idx,,drop = F]
  other_direction <- ifelse(direction == "<=", ">=", "<=")
  objective_in <- grad_vec
  constr_mat <- other_mat
  k <- nrow(constr_mat)

  if(direction == "<="){
    constr_ub <- rep(-tol, k)
    if(all(is.na(other_bound))) constr_lb <- rep(-Inf, k) else constr_lb <- rep(other_bound, k)
  } else {
    constr_lb <- rep(tol, k)
    if(all(is.na(other_bound))) constr_ub <- rep(Inf, k) else constr_ub <- rep(other_bound, k)
  }

  if(all(is.na(other_bound))){
    var_ub <- rep(Inf, k); var_lb <- rep(-Inf, k)
  } else {
    var_ub <- rep(abs(other_bound), k); var_lb <- rep(-abs(other_bound), k)
  }

  res <- clplite::clp_solve(objective_in, constr_mat, constr_lb, constr_ub, var_lb, var_ub, max = F)

  stopifnot(res$status == 0)

  if(direction == "<="){
    stopifnot(all(other_mat %*% res$solution <= 0))
  } else {
    stopifnot(all(other_mat %*% res$solution >= 0))
  }

  res$solution
}
