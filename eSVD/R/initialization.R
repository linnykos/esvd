#' Initialization for matrix factorization
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param k positive integer
#' @param family either \code{"gaussian"} or \code{"exponential"}
#' @param max_val maximum value of the inner product (with the correct sign)
#' @param max_iter numeric
#' @param tol numeric
#' @param verbose boolean
#' @param ... extra arguments
#'
#' @return list
#' @export
initialization <- function(dat, k = 2, family,
                           max_val = NA,
                           max_iter = 10, tol = 1e-3,
                           verbose = F, ...){
  direction <- .dictate_direction(family)
  if(!is.na(max_val)){
    if(!is.na(direction) && direction == "<=") max_val <- -max_val
  }

  # initialize
  dat <- .matrix_completion(dat, k = k)
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

  # projected gradient descent
  nat_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
                                          direction = direction,
                                          max_iter = max_iter,
                                          tol = tol, ...)

  res <- .svd_projection(nat_mat, k = k, factors = T)
  u_mat <- res$u_mat; v_mat <- res$v_mat

  if(!is.na(direction)){
    if(direction == "<=") {
      stopifnot(all(nat_mat[which(!is.na(dat))] < 0))
    } else {
      stopifnot(all(nat_mat[which(!is.na(dat))] > 0))
    }
  }

  tmp <- .fix_rank_defficiency_initialization(u_mat, v_mat, direction)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

  tmp <- .reparameterize(u_mat, v_mat)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

  list(u_mat = u_mat, v_mat = v_mat)
}

##################################

# enforces all the resulting entries to be non-negative
.matrix_completion <- function(dat, k){
  if(any(is.na(dat))){
    lambda0_val <- softImpute::lambda0(dat)
    res <- softImpute::softImpute(dat, rank.max = k, lambda = min(30, lambda0_val/100))
    if(k == 1){
      diag_mat <- matrix(res$d[1], 1, 1)
    } else {
      diag_mat <- diag(res$d)
    }
    pred_naive <- res$u %*% diag_mat %*% t(res$v)
    dat[which(is.na(dat))] <- pred_naive[which(is.na(dat))]
  }

  pmax(dat, 0)
}

.determine_initial_matrix <- function(dat, family, k, max_val = NA, tol = 1e-3, ...){
  dat[which(dat <= tol)] <- tol/2
  nat_mat <- .mean_transformation(dat, family, ...)
  direction <- .dictate_direction(family)

  .project_rank_feasibility(nat_mat, k = k, direction = direction,
                                    max_val = max_val)$matrix
}

.projected_gradient_descent <- function(dat, k = 2,
                                        max_val = NA, direction = "<=",
                                        max_iter = 50, tol = 1e-3,
                                        ...){
  nat_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val, ...)
  iter <- 1
  new_obj <- .evaluate_objective_mat(dat, nat_mat, ...)
  old_obj <- Inf

  while(abs(new_obj - old_obj) > tol & iter < max_iter){
    old_obj <- new_obj
    gradient_mat <- .gradient_mat(dat, nat_mat, ...)
    new_mat <- .adaptive_gradient_step(dat, nat_mat, gradient_mat, k = k,
                                       max_val = max_val, direction = direction,
                                       ...)

    new_obj <- .evaluate_objective_mat(dat, new_mat, ...)
    nat_mat <- new_mat
    iter <- iter + 1
  }

  nat_mat
}

.svd_projection <- function(mat, k, factors = F,
                            u_alone = F, v_alone = F){
  stopifnot(min(dim(mat)) >= k)

  if(min(dim(mat)) > k+2){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      RSpectra::svds(mat, k = k + 2)
    }, error = function(e){
      svd(mat)
    })
  } else {
    res <- svd(mat)
  }


  if(k == 1){
    diag_mat <- matrix(res$d[1], 1, 1)
  } else {
    diag_mat <- diag(res$d[1:k])
  }

  if(factors){
    if(u_alone){
      list(u_mat = res$u[,1:k,drop = F],
           v_mat = res$v[,1:k,drop = F]%*%diag_mat)
    } else if(v_alone) {
      list(u_mat = res$u[,1:k,drop = F]%*%diag_mat,
           v_mat = res$v[,1:k,drop = F])
    } else {
      list(u_mat = res$u[,1:k,drop = F]%*%sqrt(diag_mat),
           v_mat = res$v[,1:k,drop = F]%*%sqrt(diag_mat))
    }
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
#' @param nat_mat \code{n} by \code{d} matrix
#' @param gradient_mat \code{n} by \code{d} matrix
#' @param k numeric
#' @param max_val numeric or \code{NA}
#' @param direction "<=" or ">=" or NA
#' @param stepsize_init numeric
#' @param stepdown_factor numeric
#' @param max_iter numeric
#' @param ... other parameters
#'
#' @return \code{n} by \code{d} matrix
.adaptive_gradient_step <- function(dat, nat_mat, gradient_mat, k,
                                    max_val = NA, direction = "<=",
                                    stepsize_init = 100, stepdown_factor = 2,
                                    max_iter = 20, ...){
  stepsize <- stepsize_init
  init_obj <- .evaluate_objective_mat(dat, nat_mat, ...)
  iter <- 1

  while(iter > max_iter){
    res <- nat_mat - stepsize*gradient_mat
    new_mat <- .project_rank_feasibility(res, direction = direction,
                               max_val = max_val)$matrix

    if(!any(is.na(new_mat))){
      new_obj <- .evaluate_objective_mat(dat, new_mat, ...)

      if(new_obj < init_obj) return(new_mat)
    }

    stepsize <- stepsize/stepdown_factor
    iter <- iter + 1
  }

  # was not able to project
  nat_mat
}

# alternating projection heuristic to find intersection of two sets
.project_rank_feasibility <- function(mat, k, direction, max_val = NA,
                                      max_iter = 20, tol = 1e-6){
  stopifnot(!is.na(max_val) | !is.na(direction))
  if(!is.na(max_val) & !is.na(direction)) stopifnot((direction == "<=" & max_val < 0) | (direction == ">=" & max_val > 0))

  mat_org <- mat
  iter <- 1

  while(iter < max_iter){
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

  # if the alternating projection strategy above failed, use a SBM-projection
  mat <- .absolute_threshold(mat_org, direction, max_val)
  res <- .sbm_projection(mat, k)

  list(matrix = res, iter = NA)
}

.absolute_threshold <- function(mat, direction, max_val = NA, tol2 = 1e-3){
  if(is.na(direction)){
    max_val <- abs(max_val)

    idx <- which(abs(mat) >= max_val)
    val <- mat[idx]
    mat[idx] <- sign(val)*(max_val-tol2)

  } else if (direction == "<=") {

    tol <- -tol2
    if(any(mat < 0)) tol <- min(max(mat[mat < 0]), tol)
    stopifnot(tol < 0)
    mat[mat > -tol] <- tol
    if(!is.na(max_val)) mat[mat < max_val] <- max_val+tol2

  } else {

    tol <- tol2
    if(any(mat > 0)) tol <- max(min(mat[mat > 0]), tol)
    stopifnot(tol > 0)
    mat[mat < tol] <- tol
    if(!is.na(max_val)) mat[mat > max_val] <- max_val-tol2
  }

  mat
}

.fix_rank_defficiency_initialization <- function(u_mat, v_mat, direction){
  k <- ncol(u_mat)
  nat_mat <- u_mat %*% t(v_mat)
  k2 <- as.numeric(Matrix::rankMatrix(nat_mat))

  if(k != k2){
    stopifnot(k2 < k)
    sign_val <- ifelse(direction == ">=", 1, -1)

    sd_val <- mean(c(apply(u_mat[,1:k2, drop = F], 2, stats::sd),apply(v_mat[,1:k2, drop = F], 2, stats::sd)))
    for(i in (k2+1):k){
      u_mat[,i] <- abs(stats::rnorm(nrow(u_mat), sd = sd_val))
      v_mat[,i] <- sign_val*abs(stats::rnorm(nrow(v_mat), sd = sd_val))
    }
  }

  list(u_mat = u_mat, v_mat = v_mat)
}
