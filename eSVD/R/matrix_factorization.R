#' Fit the factorization
#'
#' Perform alternating minimization to estimate the low-rank matrix of natural
#' parameters corresponding the observations in \code{dat}. This low-rank matrix
#' is factorized by \code{u_mat} and \code{v_mat}, of which the user needs to input
#' an initial estimate of such matrices. These initial estimates can come from
#' \code{eSVD::initialization}, but can come from any other method as long
#' as \code{ncol(u_mat)==ncol(v_mat)}, \code{nrow(u_mat)==nrow(dat)} and
#' \code{nrow(v_mat)==ncol(dat)}, and the inner products between \code{u_mat} and
#' \code{v_mat} respect the domain of the natural parameters for \code{family}.
#'
#' For more information on how the optimization works, see the documentation
#' of \code{eSVD:::.optimize_mat}, which (if you look at the lower-level functions)
#' uses Frank-Wolfe to do the inner optimizations via \code{eSVD:::.frank_wolfe}
#' and choses step-sizes via \code{eSVD:::.binary_search}.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{p} columns represent genes
#' @param u_mat initial factorization, of size \code{n} by \code{k}
#' @param v_mat initial factorization, of size \code{p} by \code{k}
#' @param max_val maximum magnitude of the inner product
#' @param family character (\code{"gaussian"}, \code{"exponential"}, \code{"poisson"}, \code{"neg_binom"},
#' or \code{"curved gaussian"})
#' @param tol small positive number to dictate the convergence of the objective function
#' @param max_iter maximum number of iterations for the algorithm
#' @param verbose boolean
#' @param return_path boolean
#' @param ncores positive integer
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return list
#' @export
fit_factorization <- function(dat, u_mat, v_mat, max_val = NA,
                               family,
                               tol = 1e-3, max_iter = 100,
                               verbose = F, return_path = F,
                               ncores = NA, ...){
  if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)
  stopifnot(length(which(dat > 0)) > 0)
  stopifnot(length(which(dat < 0)) == 0)
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])
  if(!is.na(max_val)) stopifnot(max_val >= 0)

  tmp <- .check_rank(u_mat, v_mat)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat
  k <- ncol(u_mat); n <- nrow(u_mat); p <- nrow(v_mat)

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

    nat_mat <- u_mat%*%t(v_mat)
    if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(nat_mat <= 0))
    } else { stopifnot(all(nat_mat >= 0)) }
    }

    # reparameterize
    tmp <- .reparameterize(u_mat, v_mat, equal_covariance = F)
    u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

    u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, parallelized = !is.na(ncores), ...)

    nat_mat <- u_mat%*%t(v_mat)
    if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(nat_mat <= 0))
      } else { stopifnot(all(nat_mat >= 0)) }
    }

    # reparameterize
    tmp <- .reparameterize(u_mat, v_mat, equal_covariance = F)
    u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

    v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, parallelized = !is.na(ncores), ...)

    next_obj <- .evaluate_objective(dat, u_mat, v_mat, ...)


    if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))
    if(return_path) res_list[[length(res_list)+1]] <- list(u_mat = u_mat, v_mat = v_mat)

    obj_vec <- c(obj_vec, next_obj)
  }

  tmp <- .reparameterize(u_mat, v_mat, equal_covariance = T)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

  nat_mat <- u_mat%*%t(v_mat)
  if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(nat_mat <= 0))
  } else { stopifnot(all(nat_mat >= 0)) }
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

#' Optimize a matrix, row-by-row
#'
#' Given a data matrix \code{dat} where the first element of \code{class(dat)}
#' is the associated \code{family} (i.e., tells us which negative
#' log-likelihood function we want to optimize over), optimize over each row of
#' \code{current_mat} holding the values in \code{other_mat} fixed. This function
#' determines which matrix is aligned with the rows or columns based on the setting of
#' \code{left}. If \code{left=TRUE}, then the rows of \code{current_mat}
#' correspond to the rows of \code{dat}, and if \code{left=FALSE}, then the
#' rows of \code{current_mat} correspond to the columns of \code{dat}.
#'
#' This function does a gradient descent for every row in \code{current_mat}
#' via the function \code{eSVD:::.optimize_row}. See the documentation for that function for
#' more details.
#'
#' If \code{parallelized=TRUE}, it is assumed the cores is already registered
#' via \code{doMC::registerDoMC}.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{p} columns represent genes
#' @param current_mat a matrix with \code{k} columns, and number of rows equal to either
#' \code{n} or \code{p}
#' @param other_mat a matrix with \code{k} columns, and number of rows equal to either
#' \code{n} or \code{p} (whichever one isn't what \code{current_mat} has)
#' @param left boolean
#' @param max_val maximum magnitude of the inner product
#' @param parallelized boolean
#' @param verbose boolean
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return matrix of the same dimension as \code{current_mat}
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

#' Optimize a vector
#'
#' Given a data vector \code{dat_vec} where the first element of \code{class(dat_vec)}
#' is the associated \code{family} (i.e., tells us which negative
#' log-likelihood function we want to optimize over), optimize the values in
#' \code{current_vec} holding the values in \code{other_mat} fixed.
#' This function does a gradient descent, so the values in \code{current_vec}
#' are important since they determine the "initial values" of the optimization
#' procedure as well as ensure that these initial values are feasible for the optimization
#' problem.
#'
#' This function uses Frank-Wolfe, mainly using \code{eSVD:::.frank_wolfe}.
#' This function is called by \code{eSVD:::.optimize_mat}.
#'
#' @param dat_vec a vector
#' @param current_vec a vector
#' @param other_mat a matrix with the number of rows equal to \code{length(dat_vec)} and the
#' number of columns equal to \code{length(current_vec)}
#' @param n the number of rows of \code{dat} when called on in \code{eSVD:::.optimize_mat}
#' @param p the number of columns of \code{dat} when called on in \code{eSVD:::.optimize_mat}
#' @param max_iter maximum number of iterations for the algorithm
#' @param max_val maximum magnitude of the inner product
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return a vector of length \code{current_vec}
.optimize_row <- function(dat_vec, current_vec, other_mat, n, p, max_iter = 100,
                          max_val = NA, ...){
  stopifnot(length(which(!is.na(dat_vec))) > 0, length(dat_vec) == nrow(other_mat))

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

#' Binary search for the appropriate step size
#'
#' This is a variant of the line search variant of Frank-Wolfe, where we choose
#' the step size that minimizes the decrease in the objective function.
#' In this variant, we assume that the minimizing step size lies between 0 and 1,
#' and we search for the optimal step size via binary search. We first run
#' 5 iterations of binary search. Then, in the next stage,
#' we keep running binary search if (within our search space strategy), we
#' can still find a new possible value that decreases the objective function.
#'
#' @param dat_vec a vector
#' @param current_vec a vector
#' @param step_vec a vector of the same length of \code{current_vec} (typically
#' the vector output by \code{eSVD:::.frank_wolfe})
#' @param other_mat a matrix with the number of rows equal to \code{length(dat_vec)} and the
#' number of columns equal to \code{length(current_vec)}
#' @param n the number of rows of \code{dat} when called on in \code{eSVD:::.optimize_mat}
#' @param p the number of columns of \code{dat} when called on in \code{eSVD:::.optimize_mat}
#' @param max_iter maximum number of iterations for the algorithm
#' @param ...  extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return a scalar
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
#' Solves the linear program (no formatting of text).
#' Data: g_i (vector), i from 1 to k. V (matrix), with d rows and k columns.
#' Variables: y_+1 to y_+k, y_-1 to y_-k
#' Optimization problem: min g_1*(y_+1-y_-1) + ... + g_k(y_+k-y_-k)
#' such that: V %*% (y_+(1:k) - y_-(1:k)) <= tol elementwise (for entire vector of length d)
#'            V %*% (y_+(1:k) - y_-(1:k)) >= other_bound
#'            y_+i, y_-i >= 0 for i from 1 to k
#'
#' The character \code{direction} should be related to \code{other_bound}.
#' Specifically, if \code{direction="<="}, then \code{other_bound} should be a negative number,
#' meaning all the inner products between the solution of the optimization problem and the
#' rows of \code{other_mat} should be negative.
#' Otherwise, if  \code{direction=">="}, then \code{other_bound} should be a positive number,
#' meaning all the inner products between the solution of the optimization problem and the
#' rows of \code{other_mat} should be positive.
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

#' Check the rank of the inner product between the rows of two matrices
#'
#' If the rank of \code{u_mat \%*\% t(v_mat)} is not equal to
#' \code{ncol(u_mat)} (which is equal to \code{ncol(v_mat)}), then
#' drop the last few columns of \code{u_mat} and \code{v_mat}.
#' This function assumes the source of lack-of-rank is from these dropped columns.
#'
#' @param u_mat a matrix
#' @param v_mat a matrix with the same number of columns as \code{u_mat}
#'
#' @return a list containing \code{u_mat} and \code{v_mat}
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
