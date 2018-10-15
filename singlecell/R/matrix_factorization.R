.fit_factorization <- function(dat, u_mat, v_mat, tol = 1e-3,
                               max_val = NA,
                               max_iter = 100, verbose = F,
                               family = "exponential",
                               cores = NA){
  if(!is.na(cores)) doMC::registerDoMC(cores = cores)
  stopifnot(length(which(dat > 0)) > 0)
  stopifnot(length(which(dat < 0)) == 0)
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))
  k <- ncol(u_mat)
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

  idx <- which(dat == 0)
  min_val <- min(dat[which(dat > 0)])
  dat[which(dat == 0)] <- min_val/2

  current_obj <- Inf
  next_obj <- .evaluate_objective(dat, u_mat, v_mat)
  obj_vec <- c(next_obj)
  if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))

  while((is.na(tol) | abs(current_obj - next_obj) > tol) & length(obj_vec) < max_iter){
    current_obj <- next_obj

    u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, !is.na(cores))
    v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, !is.na(cores))

    next_obj <- .evaluate_objective(dat, u_mat, v_mat)

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

.evaluate_objective <- function (dat, u_mat, v_mat) {
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))

  UseMethod(".evaluate_objective")
}

.evaluate_objective.default <- function(dat, u_mat, v_mat){
  .evaluate_objective.exponential(dat, u_mat, v_mat)
}

.evaluate_objective_single <- function (dat_vec, current_vec, other_mat) {
  stopifnot(!is.matrix(dat_vec))
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  UseMethod(".evaluate_objective_single")
}

.evaluate_objective_single.default <- function(dat_vec, current_vec, other_mat){
  .evaluate_objective_single.exponential(dat_vec, current_vec, other_mat)
}

.gradient_vec <- function (dat_vec, current_vec, other_mat) {
  stopifnot(!is.matrix(dat_vec))
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  UseMethod(".gradient_vec")
}

.gradient_vec.default <- function(dat_vec, current_vec, other_mat){
  .gradient_vec.exponential(dat_vec, current_vec, other_mat)
}

#########

.optimize_mat <- function(dat, current_mat, other_mat, left = T, max_val = NA,
                          parallelized = F){
  stopifnot(length(class(dat)) == 2)

  stopifnot(ncol(current_mat) == ncol(other_mat))
  if(left) {
    stopifnot(ncol(dat) == nrow(other_mat))
  } else {
    stopifnot(nrow(dat) == nrow(other_mat))
  }

  if(parallelized){
    func <- function(i){
      if(left) { dat_vec <- dat[i,] } else { dat_vec <- dat[,i] }
      class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
      .optimize_row(dat_vec, current_mat[i,], other_mat, max_val = max_val)
    }

    lis <- foreach::"%dopar%"(foreach::foreach(i = 1:nrow(current_mat)), func(i))
    current_mat <- do.call(rbind, lis)

  } else {
    for(i in 1:nrow(current_mat)){
      if(left) { dat_vec <- dat[i,] } else { dat_vec <- dat[,i] }
      class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
      if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                                other_mat, max_val = max_val)
    }
  }


  current_mat
}

.optimize_row <- function(dat_vec, current_vec, other_mat, max_iter = 100,
                          max_val = NA){
  stopifnot(length(which(!is.na(dat_vec))) > 0)

  if(class(dat_vec)[1] == "exponential"){
    direction = "<="
  } else if(class(dat_vec)[1] == "gaussian"){
    direction = ">="
  } else {
    stop("input vector in .optimize_row() does not have a proper class")
  }

  current_obj <- Inf
  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)
  iter <- 1

  while(abs(current_obj - next_obj) > 1e-6 & iter < max_iter){
    current_obj <- next_obj

    grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat)
    step_vec <- .frank_wolfe(grad_vec, other_mat, which(!is.na(dat_vec)),
                             direction = direction, other_bound = max_val)
    step_size <- .binary_search(dat_vec, current_vec, step_vec, other_mat)
    current_vec <- (1-step_size)*current_vec + step_size*step_vec

    next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)

    iter <- iter + 1
  }

  current_vec
}

.binary_search <- function(dat_vec, current_vec, step_vec, other_mat,
                           max_iter = 100){
  form_current <- function(s){
    stopifnot(0<=s, s<=1)
    (1-s)*current_vec + s*step_vec
  }

  upper <- 1; lower <- 0; mid <- 0.5
  iter <- 1

  upper_val <- .evaluate_objective_single(dat_vec, form_current(upper), other_mat)
  lower_val <- .evaluate_objective_single(dat_vec, form_current(lower), other_mat)

  upper_val_org <- upper_val; lower_val_org <- lower_val
  mid_val <- 2*max(upper_val, lower_val);

  # stage 1: do at least a few iterations first
  while(iter <= 6){
    mid_val <- .evaluate_objective_single(dat_vec, form_current(mid), other_mat)

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
    mid_val <- .evaluate_objective_single(dat_vec, form_current(mid), other_mat)

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
  stopifnot(all(other_mat %*% res$solution <= 0))

  res$solution
}

#' L1 projection of the product of two vectors to become nonpositive
#'
#' Solves the linear program:
#' Data: u_i (vector), i from 1 to k. V (matrix), with d rows and k columns.
#' Variables: s_1 to s_k, u_+1 to u_+k, and u_-1 to u_-k.
#' Optimization problem: min s_1 + ... + s_k
#' such that: s_i >= u_i - (u_+i - u_-i) for i from 1 to k,
#'            s_i >= -(u_i - (u_+i - u_-i)) for i from 1 to k,
#'            V %*% (u_+(1:k) - u_-(1:k)) <= tol elementwise (for entire vector of length d)
#'            u_+i, u_-i, s_i >= 0 for i from 1 to k
#'
#' @param current_vec vector
#' @param other_mat matrix
#' @param idx row indices for other_mat
#' @param tol numeric
#' @param direction character
#' @param other_bound numeric
#'
#' @return
.projection_l1 <- function(current_vec, other_mat, idx = 1:nrow(other_mat),
                           tol = 0.001, direction = "<=", other_bound = NA){
  if(length(idx) == 0) return(current_vec)
  stopifnot(ncol(other_mat) == length(current_vec))
  stopifnot(is.na(other_bound) || (direction == "<=" & other_bound < 0) ||
              (direction == ">=" & other_bound > 0))

  other_mat <- other_mat[idx,,drop = F]
  if(direction == "<=") {
    if(all(other_mat %*% current_vec < 0) & (is.na(other_bound) || all(other_mat %*% current_vec > other_bound))) return(current_vec)
    other_direction <- ">="
  } else {
    if(all(other_mat %*% current_vec > 0) & (is.na(other_bound) || all(other_mat %*% current_vec < other_bound))) return(current_vec)
    other_direction <- "<="
  }

  k <- length(current_vec)
  d <- nrow(other_mat)

  objective_in <- c(rep(0, 2*k), rep(1, k))
  const_mat <- rbind(cbind(diag(k), -diag(k), diag(k)),
                     cbind(-diag(k), diag(k), diag(k)))
  const_mat2 <- t(apply(other_mat, 1, function(v){
    c(v, -v, rep(0,k))
  }))
  const_mat <- rbind(const_mat, const_mat2)
  if(!is.na(other_bound)) const_mat <- rbind(const_mat, const_mat2)

  const_dir <- c(rep(">=", 2*k), rep(direction, nrow(other_mat)))
  if(!is.na(other_bound)) const_dir <- c(const_dir, rep(other_direction, nrow(other_mat)))

  if(direction == "<=") s <- -1 else s <- 1
  const_rhs <- c(current_vec, -current_vec, rep(s*tol, nrow(other_mat)))
  if(!is.na(other_bound)) const_rhs <- c(const_rhs, rep(other_bound, nrow(other_mat)))

  res <- lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)
  res$solution[1:k] - res$solution[(k+1):(2*k)]
}
