.fit_exponential_factorization <- function(dat, k = 2, tol = 1e-3,
                                           max_iter = 100, verbose = F){
  stopifnot(length(which(dat > 0)) > 0)
  stopifnot(length(which(dat < 0)) == 0)

  idx <- which(dat == 0)
  min_val <- min(dat[which(dat > 0)])
  dat[which(dat == 0)] <- min_val/2

  init <- .initialization(dat, k = k)
  u_mat <- init$u_mat; v_mat <- init$v_mat

  current_obj <- Inf
  next_obj <- .evaluate_objective(dat, u_mat, v_mat)
  iter <- 1
  if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))

  while(current_obj - next_obj > tol & iter < max_iter){
    current_obj <- next_obj

    u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T)
    v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F)

    next_obj <- .evaluate_objective(dat, u_mat, v_mat)

    if(verbose) print(paste0("Iter ", iter, ": Decrease is ", abs(current_obj - next_obj)))

    iter <- iter + 1
  }

  pred_mat <- u_mat %*% t(v_mat)
  svd_res <- svd(pred_mat)

  if(k == 1) diag_vec <- as.matrix(sqrt(svd_res$d[1:k])) else diag_vec <- sqrt(svd_res$d[1:k])
  list(u_mat = svd_res$u[,1:k,drop = F] %*% diag(diag_vec),
       v_mat = svd_res$v[,1:k,drop = F] %*% diag(diag_vec),
       iter = iter - 1)
}

.evaluate_objective <- function(dat, u_mat, v_mat){
  stopifnot(nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat[idx] < 0))

  sum(-log(-pred_mat[idx]) - pred_mat[idx]*dat[idx])
}

.evaluate_objective_single <- function(dat_vec, current_vec, other_mat){
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] < 0))

  sum(-log(-pred_vec[idx]) - pred_vec[idx]*dat_vec[idx])
}

########


#' Initialization function
#'
#' @param dat matrix
#' @param k numeric
#'
#' @return list
#'
#' @importClassesFrom recommenderlab realRatingMatrix
.initialization <- function(dat, k = 2){
  stopifnot(length(which(dat == 0)) == 0)
  dat2 <- dat
  dat2[which(dat > 0)] <- -1/dat2[which(dat > 0)]

  if(any(is.na(dat))){
    dat2 <- methods::new("realRatingMatrix", data = recommenderlab::dropNA(dat2))

    funk_svd <- recommenderlab::funkSVD(dat2, k = k)
    u_mat <- funk_svd$U
    v_mat <- funk_svd$V

    prod_mat <- u_mat %*% t(v_mat)
    svd_res <- svd(prod_mat)

    if(k == 1) diag_vec <- as.matrix(sqrt(svd_res$d[1:k])) else diag_vec <- sqrt(svd_res$d[1:k])
    u_mat <- svd_res$u[,1:k,drop = F] %*% diag(diag_vec)
    v_mat <- svd_res$v[,1:k,drop = F] %*% diag(diag_vec)
  } else {
    svd_res <- svd(dat2)

    if(k == 1) diag_vec <- as.matrix(sqrt(svd_res$d[1:k])) else diag_vec <- sqrt(svd_res$d[1:k])
    u_mat <- svd_res$u[,1:k,drop = F] %*% diag(diag_vec)
    v_mat <- svd_res$v[,1:k,drop = F] %*% diag(diag_vec)
  }

  # project v back into positive space based on u
  for(j in 1:nrow(v_mat)){
    v_mat[j,] <- .projection_l1(v_mat[j,], u_mat, which(!is.na(dat[,j])))
  }

  list(u_mat = u_mat, v_mat = v_mat)
}

#########

.optimize_mat <- function(dat, current_mat, other_mat, left = T){
  stopifnot(ncol(current_mat) == ncol(other_mat))
  if(left) {
    stopifnot(ncol(dat) == nrow(other_mat))
  } else {
    stopifnot(nrow(dat) == nrow(other_mat))
  }

  for(i in 1:nrow(current_mat)){
    if(left) {dat_vec <- dat[i,]} else {dat_vec <- dat[,i]}
    if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,], other_mat)
  }

  current_mat
}

.optimize_row <- function(dat_vec, current_vec, other_mat, max_iter = 100){
  stopifnot(length(which(!is.na(dat_vec))) > 0)

  current_obj <- Inf
  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)
  iter <- 1

  while(abs(current_obj - next_obj) > 1e-6 & iter < max_iter){
    current_obj <- next_obj

    grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat)
    stepsize <- .backtrack_linesearch(dat_vec, current_vec, other_mat, grad_vec)

    current_vec <- current_vec - stepsize * grad_vec
    current_vec <- .projection_l1(current_vec, other_mat, which(!is.na(dat_vec)))

    next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)

    iter <- iter + 1
  }

  current_vec
}

.gradient_vec <- function(dat_vec, current_vec, other_mat){
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(-1/pred_vec[j]-dat_vec[j])
  })

  if(is.matrix(tmp)) rowSums(tmp) else sum(tmp)
}

.backtrack_linesearch <- function(dat_vec, current_vec, other_mat, grad_vec,
                                  t_init = 10, beta = .5, alpha = .5){
  t_current <- t_init
  idx <- which(!is.na(dat_vec))

  while(TRUE){
    if(any((other_mat %*% (current_vec - t_current*grad_vec))[idx] >= 0)){
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
#'            V %*% (u_+(1:k) - u_-(1:k)) <= 0 elementwise (for entire vector of length d)
#'            s_i, u_+i, u_-i >= 0 for i from 1 to k
#'
#' @param current_vec vector
#' @param other_mat matrix
#' @param idx row indices for other_mat
#' @param tol numeric
#'
#' @return
.projection_l1 <- function(current_vec, other_mat, idx = 1:nrow(other_mat),
                           tol = 1e-6){
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

  const_dir <- c(rep(">=", 2*k), rep("<=", nrow(other_mat)))
  const_rhs <- c(current_vec, -current_vec, rep(-tol, nrow(other_mat)))

  res <- lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)
  res$solution[1:k] - res$solution[(k+1):(2*k)]
}
