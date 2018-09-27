.fit_gaussian_factorization <- function(dat, k = 2, lambda = 0.01,
                                           max_iter = 100){
  init <- .initialization(dat)
  u_mat <- init$u_mat; v_mat <- init$v_mat

  current_obj <- Inf
  next_obj <- .evaluate_objective(dat, u_mat, v_mat)
  iter <- 1

  while(abs(current_obj - next_obj) > 1e-6 & iter < max_iter){
    current_obj <- next_obj

    u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T)
    v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F)

    next_obj <- .evaluate_objective(dat, u_mat, v_mat)

    iter <- iter + 1
  }

  pred_mat <- u_mat %*% t(v_mat)
  svd_res <- svd(pred_mat)

  list(u_mat = svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k])),
       v_mat = svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k])),
       iter = iter - 1)
}

.evaluate_objective <- function(dat, u_mat, v_mat){
  stopifnot(nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))

  sum(pred_mat[idx]*(1-dat[idx]))
}

.evaluate_objective_single <- function(dat_vec, current_vec, other_mat){
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  sum(pred_vec[idx]*(1-dat_vec[idx]))
}

########


#' Initialization function
#'
#' @param dat matrix
#' @param k numeric
#' @param lambda numeric
#'
#' @return list
#'
#' @importClassesFrom recommenderlab realRatingMatrix
.initialization <- function(dat, k = 2, lambda = 0.01){
  dat2 <- methods::new("realRatingMatrix", data = recommenderlab::dropNA(dat))

  funk_svd <- recommenderlab::funkSVD(dat2, k = k, lambda = lambda)
  u_mat <- funk_svd$U
  v_mat <- funk_svd$V

  prod_mat <- u_mat %*% t(v_mat)
  svd_res <- svd(prod_mat)

  u_mat <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
  v_mat <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

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
    next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)

    iter <- iter + 1
  }

  current_vec
}

.gradient_vec <- function(dat_vec, current_vec, other_mat){
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  non_na_idx <- which(!is.na(dat_vec))
  tmp <- sapply(non_na_idx, function(j){
    other_mat[j,]*(1-dat_vec[j])
  })

  rowSums(tmp)
}

.backtrack_linesearch <- function(dat_vec, current_vec, other_mat, grad_vec,
                                  beta = .5, alpha = .5, tol = 1e-4){
  t_current <- min(as.numeric(other_mat %*% current_vec)/as.numeric(other_mat %*% grad_vec))
  t_current <- max(t_current - tol, t_current/2)

  while(TRUE){
    obj1 <- .evaluate_objective_single(dat_vec, current_vec - t_current*grad_vec, other_mat)
    obj2 <- .evaluate_objective_single(dat_vec, current_vec, other_mat) - alpha*t_current*.l2norm(grad_vec)^2

    if(obj1 > obj2) t_current <- t_current*beta else break()
  }

  t_current
}

.l2norm <- function(x){sqrt(sum(x^2))}

#######################

.projection_l1 <- function(current_vec, other_mat, idx = 1:nrow(other_mat)){
  if(length(idx) == 0) return(current_vec)
  stopifnot(ncol(other_mat) == length(current_vec))

  other_mat <- other_mat[idx,]

  k <- length(current_vec)
  d <- nrow(other_mat)

  objective_in <- c(rep(0, 2*k), rep(1, k))
  const_mat <- rbind(cbind(diag(k), -diag(k), diag(k)),
                     cbind(-diag(k), diag(k), diag(k)))
  const_mat2 <- t(apply(other_mat, 1, function(v){
    c(v, -v, rep(0,k))
  }))
  const_mat <- rbind(const_mat, const_mat2)

  const_dir <- c(rep(">=", 2*k), rep(">=", nrow(other_mat)))
  const_rhs <- c(current_vec, -current_vec, rep(0, nrow(other_mat)))

  res <- lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)
  res$solution[1:k] - res$solution[(k+1):(2*k)]
}
