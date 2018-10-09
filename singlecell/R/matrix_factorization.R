.fit_factorization <- function(dat, u_mat, v_mat, tol = 1e-3,
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

  while(current_obj - next_obj > tol & length(obj_vec) < max_iter){
    current_obj <- next_obj

    u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, !is.na(cores))
    v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, !is.na(cores))

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

.optimize_mat <- function(dat, current_mat, other_mat, left = T, parallelized = F){
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
      .optimize_row(dat_vec, current_mat[i,], other_mat)
    }

    lis <- foreach::"%dopar%"(foreach::foreach(i = 1:nrow(current_mat)), func(i))
    current_mat <- do.call(rbind, lis)

  } else {
    for(i in 1:nrow(current_mat)){
      if(left) { dat_vec <- dat[i,] } else { dat_vec <- dat[,i] }
      class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
      if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,], other_mat)
    }
  }


  current_mat
}

.optimize_row <- function(dat_vec, current_vec, other_mat, max_iter = 100){
  stopifnot(length(which(!is.na(dat_vec))) > 0)

  if(class(dat_vec)[1] == "exponential"){
    direction = "<="
  } else if(class(dat_vec)[1] == "gaussian"){
    direction = ">="
  }

  current_obj <- Inf
  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)
  iter <- 1

  while(abs(current_obj - next_obj) > 1e-6 & iter < max_iter){
    current_obj <- next_obj

    grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat)
    stepsize <- .backtrack_linesearch(dat_vec, current_vec, other_mat, grad_vec)

    current_vec <- current_vec - stepsize * grad_vec
    current_vec <- .projection_l1(current_vec, other_mat, which(!is.na(dat_vec)),
                                  direction = direction)

    next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)

    iter <- iter + 1
  }

  current_vec
}

.backtrack_linesearch <- function(dat_vec, current_vec, other_mat, grad_vec,
                                  t_init = 10, beta = .5, alpha = .5){
  t_current <- t_init
  idx <- which(!is.na(dat_vec))

  s <- -1*unique(as.numeric(sign(other_mat %*% current_vec)))
  stopifnot(any(s %in% c(-1,1)), !all(c(-1,1) %in% s))
  if(length(s) > 1) s <- s[s %in% c(-1,1)]

  while(TRUE){
    if(any((s*(other_mat %*% (current_vec - t_current*grad_vec))[idx]) >= 0)){
      t_current <- t_current*beta
    } else {
      obj1 <- .evaluate_objective_single(dat_vec, current_vec - t_current*grad_vec, other_mat)
      obj2 <- .evaluate_objective_single(dat_vec, current_vec, other_mat) - alpha*t_current*.l2norm(grad_vec)^2
      if(obj1 > obj2) t_current <- t_current*beta else break()
    }
  }

  t_current
}

.l2norm <- function(x){sqrt(sum(x^2))}

#######################

#' L1 projection of the product of two vectors to become nonpositive
#'
#' Solves the linear program:
#' Data: u_i (vector), i from 1 to k. V (matrix), with d rows and k columns.
#' Variables: s_1 to s_k, u_+1 to u_+k, and u_-1 to u_-k.
#' Optimization problem: min s_1 + ... + s_k
#' such that: s_i >= u_i - (u_+i - u_-i) for i from 1 to k,
#'            s_i >= -(u_i - (u_+i - u_-i)) for i from 1 to k,
#'            V %*% (u_+(1:k) - u_-(1:k)) <= tol elementwise (for entire vector of length d)
#'            s_i, u_+i, u_-i >= 0 for i from 1 to k
#'
#' @param current_vec vector
#' @param other_mat matrix
#' @param idx row indices for other_mat
#' @param tol numeric
#'
#' @return
.projection_l1 <- function(current_vec, other_mat, idx = 1:nrow(other_mat),
                           tol = 1e-6, direction = "<="){
  if(length(idx) == 0) return(current_vec)
  stopifnot(ncol(other_mat) == length(current_vec))

  other_mat <- other_mat[idx,,drop = F]
  if(all(other_mat %*% current_vec <= 0)) return(current_vec)

  k <- length(current_vec)
  d <- nrow(other_mat)

  objective_in <- c(rep(0, 2*k), rep(1, k))
  const_mat <- rbind(cbind(diag(k), -diag(k), diag(k)),
                     cbind(-diag(k), diag(k), diag(k)))
  const_mat2 <- t(apply(other_mat, 1, function(v){
    c(v, -v, rep(0,k))
  }))
  const_mat <- rbind(const_mat, const_mat2)

  const_dir <- c(rep(">=", 2*k), rep(direction, nrow(other_mat)))
  if(direction == "<=") s <- -1 else s <- 1
  const_rhs <- c(current_vec, -current_vec, rep(s*tol, nrow(other_mat)))

  res <- lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)
  res$solution[1:k] - res$solution[(k+1):(2*k)]
}
