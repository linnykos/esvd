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
                            max_val = NA, verbose = F){
  stopifnot(length(extra_weights) == nrow(dat))

  # complete the matrix
  if(any(is.na(dat))){
    lambda0_val <- softImpute::lambda0(dat)
    res <- softImpute::softImpute(dat, rank.max = k, lambda = min(30, lambda0_val/100))
    pred_naive <- res$u %*% diag(res$d) %*% t(res$v)
    dat[which(is.na(dat))] <- pred_naive[which(is.na(dat))]
  }

  idx <- which(dat == 0)
  min_val <- min(dat[which(dat > 0)])
  dat[which(dat <= 0)] <- min_val/2
  direction <- .dictate_direction(family)
  dat2 <- dat
  for(i in 1:nrow(dat2)){
    dat2[i,] <- dat2[i,]/extra_weights[i]
  }
  dat2 <- .mean_transformation(dat2, family)

  svd_res <- svd(dat2)

  if(k == 1) diag_vec <- as.matrix(sqrt(svd_res$d[1:k])) else diag_vec <- sqrt(svd_res$d[1:k])
  u_mat <- svd_res$u[,1:k,drop = F] %*% diag(diag_vec)
  v_mat <- svd_res$v[,1:k,drop = F] %*% diag(diag_vec)

  # project v back into positive space based on u
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
    u_mat <- pmax(u_mat, 1e-3)

    for(j in 1:nrow(v_mat)){
      v_mat[j,] <- .projection_l1(v_mat[j,], u_mat, direction = direction,
                                  other_bound = max_val)
    }
  }

  pred_mat <- u_mat %*% t(v_mat)
  if(direction == "<=") {
    stopifnot(all(pred_mat[which(!is.na(dat))] < 0))
  } else {
    stopifnot(all(pred_mat[which(!is.na(dat))] > 0))
  }

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

