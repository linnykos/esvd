#' Fit the factorization
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param u_mat initial factorization, of size \code{n} by \code{k}
#' @param v_mat initial factorization, of size \code{d} by \code{k}
#' @param max_val maximum value of the inner product (with the correct sign)
#' @param family character such as \code{"gaussian"} or \code{"exponential"}
#' @param tol small positive number to dictate the convergence of the objective function
#' @param max_iter maximum number of iterations for the algorithm
#' @param verbose boolean
#' @param return_path boolean
#' @param cores positive integer
#' @param ... additional parameters for the distribution
#'
#' @return list
#' @export
fit_factorization <- function(dat, u_mat, v_mat, max_val = NA,
                               family,
                               tol = 1e-3, max_iter = 100,
                               verbose = F, return_path = F,
                               cores = NA, ...){
  if(!is.na(cores)) doMC::registerDoMC(cores = cores)
  stopifnot(length(which(dat > 0)) > 0)
  stopifnot(length(which(dat < 0)) == 0)
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])
  if(!is.na(max_val)) stopifnot(max_val >= 0)

  tmp <- .check_rank(u_mat, v_mat)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat
  k <- ncol(u_mat)

  idx <- which(!is.na(dat))
  min_val <- min(dat[which(dat > 0)])
  dat[which(dat == 0)] <- min_val/2
  direction <- .dictate_direction(class(dat)[1])
  if(!is.na(direction) && direction == "<=" && !is.na(max_val)) max_val <- -max_val

  current_obj <- Inf
  next_obj <- .evaluate_objective(dat, u_mat, v_mat, ...)
  obj_vec <- c(next_obj)
  if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))
  if(return_path) res_list <- list(list(u_mat = u_mat, v_mat = v_mat)) else res_list <- NA

  while((is.na(tol) | abs(current_obj - next_obj) > tol) & length(obj_vec) < max_iter){
    current_obj <- next_obj

    pred_mat <- u_mat%*%t(v_mat)
    if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(pred_mat <= 0))
    } else { stopifnot(all(pred_mat >= 0)) }
    }

    # reparameterize
    tmp <- .reparameterize(u_mat, v_mat)
    u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

    u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, parallelized = !is.na(cores), ...)

    pred_mat <- u_mat%*%t(v_mat)
    if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(pred_mat <= 0))
      } else { stopifnot(all(pred_mat >= 0)) }
    }

    # reparameterize
    tmp <- .reparameterize(u_mat, v_mat)
    u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

    v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, parallelized = !is.na(cores), ...)

    next_obj <- .evaluate_objective(dat, u_mat, v_mat, ...)


    if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))
    if(return_path) res_list[[length(res_list)+1]] <- list(u_mat = u_mat, v_mat = v_mat)

    obj_vec <- c(obj_vec, next_obj)
  }

  tmp <- .reparameterize(u_mat, v_mat)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

  pred_mat <- u_mat%*%t(v_mat)
  if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(pred_mat <= 0))
  } else { stopifnot(all(pred_mat >= 0)) }
  }

  list(u_mat = u_mat, v_mat = v_mat, obj_vec = obj_vec, res_list = res_list)
}

#########

.diag_matrix <- function(vec){
  k <- length(vec)
  if(k == 1) {
    matrix(vec, 1, 1)
  } else {
    diag(vec)
  }
}

.optimize_mat <- function(dat, current_mat, other_mat, left = T, max_val = NA,
                          parallelized = F, verbose = F, ...){
  stopifnot(length(class(dat)) == 2)
  n <- nrow(dat); p <- ncol(dat)

  stopifnot(ncol(current_mat) == ncol(other_mat))
  if(left) {
    stopifnot(ncol(dat) == nrow(other_mat))
  } else {
    stopifnot(nrow(dat) == nrow(other_mat))
  }

  if(parallelized){
    func <- function(i){
      if(left) {
        dat_vec <- dat[i,]
      } else {
        dat_vec <- dat[,i]
      }
      class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
      .optimize_row(dat_vec, current_mat[i,], other_mat, max_val = max_val,
                    n = n, p = p, ...)
    }

    lis <- foreach::"%dopar%"(foreach::foreach(i = 1:nrow(current_mat)), func(i))
    current_mat <- do.call(rbind, lis)

  } else {
    for(i in 1:nrow(current_mat)){
      if(verbose) print(paste0("Working on index ", i, " for ", ifelse(left, "left", "right")))

      if(left) {
        dat_vec <- dat[i,]
      } else {
        dat_vec <- dat[,i]
      }
      class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
      if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                                other_mat, max_val = max_val,
                                                                n = n, p = p, ...)
    }
  }

  current_mat
}

.optimize_row <- function(dat_vec, current_vec, other_mat, n, p, max_iter = 100,
                          max_val = NA, ...){
  stopifnot(length(which(!is.na(dat_vec))) > 0)

  direction <- .dictate_direction(class(dat_vec)[1])
  current_obj <- Inf
  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat, n = n, p = p, ...)
  iter <- 1
  while(abs(current_obj - next_obj) > 1e-6 & iter < max_iter){
    current_obj <- next_obj

    grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat, n = n, p = p, ...)
    step_vec <- .frank_wolfe(grad_vec, other_mat,
                             direction = direction, other_bound = max_val)
    step_size <- .binary_search(dat_vec, current_vec, step_vec, other_mat, n = n, p = p, ...)
    current_vec <- (1-step_size)*current_vec + step_size*step_vec

    next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat, n = n, p = p, ...)

    if(!is.na(max_val)) stopifnot(all(abs(other_mat %*% current_vec) <= abs(max_val)+1e-6))
    iter <- iter + 1
  }

  current_vec
}

.binary_search <- function(dat_vec, current_vec, step_vec, other_mat,
                           n, p,
                           max_iter = 100, ...){
  form_current <- function(s){
    stopifnot(0<=s, s<=1)
    (1-s)*current_vec + s*step_vec
  }

  upper <- 1; lower <- 0; mid <- 0.5
  iter <- 1

  upper_val <- .evaluate_objective_single(dat_vec, form_current(upper), other_mat,
                                          n = n, p = p, ...)
  lower_val <- .evaluate_objective_single(dat_vec, form_current(lower), other_mat,
                                          n = n, p = p, ...)

  upper_val_org <- upper_val; lower_val_org <- lower_val
  mid_val <- 2*max(upper_val, lower_val);

  # stage 1: do at least a few iterations first
  while(iter <= 6){
    mid_val <- .evaluate_objective_single(dat_vec, form_current(mid), other_mat,
                                          n = n, p = p, ...)

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
    mid_val <- .evaluate_objective_single(dat_vec, form_current(mid), other_mat,
                                          n = n, p = p, ...)

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
#' @param tol numierc
#' @param direction character
#' @param other_bound numeric
#'
#' @return vector
.frank_wolfe <- function(grad_vec, other_mat,
                         tol = 0.001, direction = "<=", other_bound = NA){
  stopifnot(!is.na(direction) | !is.na(other_bound))
  if(!is.na(other_bound) & !is.na(direction)) stopifnot((direction == "<=" & other_bound < 0) | (direction == ">=" & other_bound > 0))

  k <- length(grad_vec)
  if(is.na(direction)){
    other_direction <- NA
  } else if(direction == "<=") {
    other_direction <- ">="
  } else {
    other_direction <- "<="
  }

  objective_in <- grad_vec
  constr_mat <- other_mat
  k <- nrow(constr_mat)

  if (is.na(direction)){
    constr_lb <- rep(-abs(other_bound), k)
    constr_ub <- rep(abs(other_bound), k)

  } else if(direction == "<="){
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

  if(!is.na(direction)){
    if(direction == "<="){
      stopifnot(all(other_mat %*% res$solution <= 0))
    } else {
      stopifnot(all(other_mat %*% res$solution >= 0))
    }
  }


  res$solution
}

.check_rank <- function(u_mat, v_mat){
  stopifnot(ncol(u_mat) == ncol(v_mat))
  k <- ncol(u_mat)
  nat_mat <- u_mat %*% t(v_mat)

  k2 <- as.numeric(Matrix::rankMatrix(nat_mat))

  if(k != k2){
    stopifnot(k2 < k)
    warning("Initial matrices given are rank defficient -- dropping ranks")
    u_mat <- u_mat[,1:k2]
    v_mat <- v_mat[,1:k2]
  }

  return(list(u_mat = u_mat, v_mat = v_mat))
}
