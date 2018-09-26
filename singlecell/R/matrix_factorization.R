.fit_exponential_factorization <- function(dat, max_iter = 100){
  init <- .initialization(dat)
  u_mat <- init$u_mat; v_mat <- init$v_mat

  current_obj <- .evaluate_objective(dat, u_mat, v_mat)
  next_obj <- Inf
  counter <- 1

  while(abs(current_obj - next_obj) < 1e-6){
    u_mat <- .update_mat(dat, u_mat, v_mat, left = T)
    v_mat <- .update_mat(dat, v_mat, u_mat, left = T)

    next_obj <- .evaluate_objective(dat, u_mat, v_mat)

    counter <- counter + 1
  }

  return(u_mat = u_mat, v_mat = v_mat)
}

.evaluate_objective <- function(dat, u_mat, v_mat){
  stopifnot(nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))

  sum(-log(pred_mat[idx]) - dat[idx] * pred_mat[idx])
}

.evaluate_objective_single <- function(dat_vec, current_vec, other_mat){
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  sum(-log(pred_vec[idx]) - dat_vec[idx] * pred_vec[idx])
}

########

.initialization <- function(dat){

}

#########

.update_mat <- function(dat, current_mat, other_mat, left = T){
  stopifnot(length(current_vec) == ncol(other_mat))
  if(left) {
    stopifnot(ncol(dat) == nrow(other_mat))
  } else {
    stopifnot(nrow(dat) == nrow(other_mat))
  }

  for(i in 1:nrow(current_mat)){
    if(left) {dat_vec <- dat[i,]} else {dat_vec <- dat[,i]}
    grad_vec <- .gradient_vec(dat_vec, current_vec[i,], other_mat)
    stepsize <- .backtrack_linesearch(dat_vec, current_vec, other_mat, grad_vec)

    current_mat[i,] <- current_mat[i,] - stepsize * grad_vec
  }

  current_mat
}

.gradient_vec <- function(dat_vec, current_vec, other_mat){
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  non_na_idx <- which(!is.na(dat_vec))
  tmp <- sapply(non_na_idx, function(j){
    -other_mat[j,]/(current_vec %*% other_mat[j,]) - dat_vec[j] * other_mat[j,]
  })

  rowSums(tmp)
}

.backtrack_linesearch <- function(dat_vec, current_vec, other_mat,
                                  grad_vec,
                                  t_init = 1, beta = .5, alpha = .5){
  t_current <- t_init

  while(TRUE){
    obj1 <- .evaluate_objective_single(current_vec - t_current*grad_vec)
    obj2 <- .evaluate_objective_single(current_vec) - alpha*t_current*.l2norm(grad_vec)^2

    if(obj1 > obj2) t_current <- t_current*beta else break()
  }

  t_current
}

.l2norm <- function(x){sqrt(sum(x^2))}
