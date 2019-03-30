#' Initialization for matrix factorization
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param k positive integer
#' @param family either \code{"gaussian"}, \code{"exponential"} or \code{"poisson"}
#' @param extra_weights vector of weights, of length \code{n}
#' @param max_val maximum value of the inner product (with the correct sign)
#' @param verbose boolean
#'
#' @return list
#' @export
initialization <- function(dat, k = 2, family = "exponential",
                           extra_weights = rep(1, nrow(dat)),
                           max_val = NA,
                           max_iter = 10, tol = 1e-3,
                           verbose = F, ...){
  stopifnot(length(extra_weights) == nrow(dat))
  direction <- .dictate_direction(family)

  # initialize
  dat <- .matrix_completion(dat, k = k)
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

  # projected gradient descent
  pred_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
                                          direction = direction,
                                          max_iter = max_iter,
                                          tol = tol, ...)

  res <- .svd_projection(pred_mat, k = k, factors = T)
  u_mat <- res$u_mat; v_mat <- res$v_mat

  if(direction == "<=") {
    stopifnot(all(pred_mat[which(!is.na(dat))] < 0))
  } else {
    stopifnot(all(pred_mat[which(!is.na(dat))] > 0))
  }

  list(u_mat = u_mat, v_mat = v_mat)
}

##################################

.matrix_completion <- function(dat, k = 2){
  if(any(is.na(dat))){
    lambda0_val <- softImpute::lambda0(dat)
    res <- softImpute::softImpute(dat, rank.max = k, lambda = min(30, lambda0_val/100))
    pred_naive <- res$u %*% diag(res$d) %*% t(res$v)
    dat[which(is.na(dat))] <- pred_naive[which(is.na(dat))]
  }

  dat
}

.projected_gradient_descent <- function(dat, k = 2,
                                        max_val = NA, direction = "<=",
                                        max_iter = 50, tol = 1e-3,
                                        ...){
  n <- nrow(dat); d <- ncol(dat)
  pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)
  iter <- 1
  new_obj <- .evaluate_objective_mat(dat, pred_mat, ...)
  old_obj <- Inf

  while(abs(new_obj - old_obj) > tol & iter < max_iter){
    old_obj <- new_obj
    gradient_mat <- .gradient_mat(dat, pred_mat, ...)
    new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = k,
                                       max_val = max_val, direction = direction,
                                       ...)

    new_obj <- .evaluate_objective_mat(dat, new_mat, ...)
    pred_mat <- new_mat
    iter <- iter + 1
  }

  pred_mat
}

.determine_initial_matrix <- function(dat, family, k, max_val = NA){
  min_val <- min(dat[which(dat > 0)])
  dat[which(dat <= 0)] <- min_val/2
  pred_mat <- .mean_transformation(dat, family)
  direction <- .dictate_direction(family)

  res <- .svd_projection(pred_mat, k = k, factors = T)
  res <- .ensure_feasibility(res$u_mat, res$v_mat, direction = direction,
                      max_val = max_val, verbose = F)
  res$u_mat %*% t(res$v_mat)
}

.svd_projection <- function(mat, k, factors = F){
  res <- svd(mat)

  if(k == 1){
    diag_mat <- matrix(res$d[1], 1, 1)
  } else {
    diag_mat <- diag(res$d[1:k])
  }

  if(factors){
    list(u_mat = res$u[,1:k,drop = F]%*%sqrt(diag_mat),
         v_mat = res$v[,1:k,drop = F]%*%sqrt(diag_mat))
  } else {
    res$u[,1:k,drop = F] %*% diag_mat %*% t(res$v[,1:k,drop = F])
  }
}

#' Adaptive projective gradient descent
#'
#' The projective gradient descent is "adaptive" in the sense
#' that it will find an appropriate step size to ensure that after projection,
#' the objective value descreases. This is a heuristic to simply enable
#' reasonable results, not necessarily theoretically justified or computationally
#' efficient.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param pred_mat \code{n} by \code{d} matrix
#' @param gradient_mat \code{n} by \code{d} matrix
#' @param k numeric
#' @param stepsize_init numeric
#' @param stepdown_factor numeric
#' @param max_iter numeric
#' @param ... other parameters
#'
#' @return \code{n} by \code{d} matrix
.adaptive_gradient_step <- function(dat, pred_mat, gradient_mat, k,
                                    max_val = NA, direction = "<=",
                                    stepsize_init = 100, stepdown_factor = 2,
                                    max_iter = 20, ...){
  stepsize <- stepsize_init
  init_obj <- .evaluate_objective_mat(dat, pred_mat, ...)
  iter <- 1

  while(TRUE){
    res <- .svd_projection(pred_mat - stepsize*gradient_mat, k = k, factors = T)
    res <- .ensure_feasibility(res$u_mat, res$v_mat, direction = direction,
                               max_val = max_val, verbose = F)
    new_mat <- res$u_mat %*% t(res$v_mat)

    new_obj <- .evaluate_objective_mat(dat, new_mat, ...)

    if(new_obj < init_obj) break() else stepsize <- stepsize/stepdown_factor
    iter <- iter + 1
    if(iter > max_iter) stop("Adaptive gradient initialization failed")
  }

  new_mat
}

.ensure_feasibility <- function(u_mat, v_mat, direction, max_val = NA,
                                verbose = F){

  for(j in 1:nrow(v_mat)){
    res <- .projection_l1(v_mat[j,], u_mat, direction = direction, other_bound = max_val)
    if(attr(res, "status") != 0) break()
    v_mat[j,] <- as.numeric(res)
  }

  if(attr(res, "status") != 0){
    if(verbose) warning("Had to use some janky fixes")
    if(length(which(u_mat < 0)) > length(which(u_mat > 0))) {
      u_mat <- -u_mat; v_mat <- -v_mat
    }
    if(direction == ">=") u_mat <- pmax(u_mat, 1e-3) else u_mat <- pmin(u_mat, -1e-3)

    for(j in 1:nrow(v_mat)){
      v_mat[j,] <- .projection_l1(v_mat[j,], u_mat, direction = direction,
                                  other_bound = max_val)
    }
  }

  if(direction == ">=") stopifnot(all(u_mat %*% t(v_mat) >= 0)) else stopifnot(all(u_mat %*% t(v_mat) <= 0))

  list(u_mat = u_mat, v_mat = v_mat)
}


#' L1 projection of the product of two vectors to become nonpositive
#'
#' Solves the linear program:
#' Data: u_i (vector), i from 1 to k. V (matrix), with d rows and k columns.
#' Variables: s_1 to s_k, z_1 to z_k
#' Optimization problem: min s_1 + ... + s_k
#' such that: s_i >= u_i - z_i for i from 1 to k,
#'            s_i >= -(u_i - z_i) for i from 1 to k,
#'            V %*% z_k <= tol elementwise (for entire vector of length d)
#'            s_i >= 0 for i from 1 to k
#'
#' @param current_vec vector
#' @param other_mat matrix
#' @param tol numeric
#' @param direction character
#' @param other_bound numeric
#'
#' @return
.projection_l1 <- function(current_vec, other_mat,
                           tol = 0.001, direction = "<=", other_bound = NA){
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

  vec <- res$solution[(k+1):(2*k)]
  attr(vec, "status") <- res$status
  vec
}

