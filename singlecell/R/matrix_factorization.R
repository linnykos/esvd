estimate_row <- function(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out, tol = 1e-5){

}

.subgradient_vec <- function(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out){
  stopifnot(ncol(latent_mat) == length(initial_vec), ncol(dat) == nrow(latent_mat))
  stopifnot(length(fixed_idx) <= nrow(dat))
  stopifnot(max(c(index_in, index_out)) <= prod(dim(dat)))
  stopifnot(length(intersect(index_in, index_out)) == 0)

  vec <- initial_vec %*% t(latent_mat[index_out,])
  vec <- sapply(vec, function(x){max(0, x)})
  second_term <- 2*vec%*% latent_mat[index_out,]

  -(dat[fixed_idx,index_in] - initial_vec %*% t(latent_mat[index_in,])) %*% latent_mat[index_in,]/length(index_in) +
    second_term/length(index_out)
}

.evaluate_objective_single <- function(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out){
  stopifnot(ncol(latent_mat) == length(initial_vec), ncol(dat) == nrow(latent_mat))
  stopifnot(length(fixed_idx) <= nrow(dat))
  stopifnot(max(c(index_in, index_out)) <= prod(dim(dat)))
  stopifnot(length(intersect(index_in, index_out)) == 0)

  vec <- initial_vec %*% t(latent_mat[index_out,])
  vec <- sapply(vec, function(x){max(0, x)})
  second_term <- sum(vec^2)

  sum((dat[fixed_idx,index_in] - initial_vec %*% t(latent_mat[index_in,]))^2)/(2*length(index_in)) +
    second_term/(2*length(index_out))
}

.evaluate_objective_full <- function(dat, u_mat, v_mat, index_in, index_out){

}
