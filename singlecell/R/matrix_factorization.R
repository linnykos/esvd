estimate_latent <- function(dat, k, dropout_func, threshold,
                            tol = 1e-5, max_iter = 500, max_outer_iter = 10, cores = 1,
                            verbose = F){
  if(verbose) print("Starting")
  res_svd <- svd(dat)
  u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
  v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

  #determine indices
  index_in_vec <- which(dat != 0)
  index_zero <- which(dat == 0)
  index_out_old <- .predict_true_zero(u_mat %*% t(v_mat), dropout_func, threshold, index_zero)
  if(length(index_out_old) == 0) stop("Threshold is too low")
  counter <- 1

  while(TRUE){
    obj_old <- .evaluate_objective_full(dat, u_mat, v_mat, index_in_vec, index_out_old)

    if(verbose) print(paste0("Starting round ", counter, " with objective value of: ", obj_old))

    while(TRUE){
      u_mat <- .estimate_matrix(dat, u_mat, v_mat, index_in_vec, index_out_old,
                                tol, max_iter, row = T, cores)
      v_mat <- .estimate_matrix(dat, v_mat, u_mat, index_in_vec, index_out_old,
                                tol, max_iter, row = F, cores)
      obj_new <- .evaluate_objective_full(dat, u_mat, v_mat, index_in_vec, index_out_old)
      if(abs(obj_old - obj_new) <= tol) break()

      if(verbose) print(obj_new)

      obj_old <- obj_new
    }

    index_out_new <- .predict_true_zero(u_mat %*% t(v_mat), dropout_func, threshold, index_zero)
    if(length(index_out_old) == length(index_out_new) && all(sort(index_out_old) == sort(index_out_new))) break()
    counter <- counter + 1
    if(counter > max_outer_iter) break()
  }

  pred_mat <- u_mat %*% t(v_mat)
  res_svd <- svd(pred_mat)

  list(u_mat = res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k])),
       v_mat = res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k])))
}

.predict_true_zero <- function(pred_dat, dropout_func, threshold, index_vec){
  prob_vec <- sapply(pred_dat[index_vec], dropout_func)
  index_vec[which(prob_vec <= threshold)]
}

.estimate_matrix <- function(dat, initial_mat, latent_mat, index_in_vec, index_out_vec,
                             tol = 1e-5, max_iter = 500, row = T, cores = 1){
  stopifnot(((row & nrow(dat) == nrow(initial_mat) & ncol(dat) == nrow(latent_mat)) |
               (!row & ncol(dat) == nrow(initial_mat) & nrow(dat) == nrow(latent_mat))))
  stopifnot(ncol(initial_mat) == ncol(latent_mat))
  stopifnot(all(index_in_vec <= prod(dim(dat))),
            all(index_out_vec <= prod(dim(dat))))

  doMC::registerDoMC(cores = cores)

  index_in_mat <- .convert_index_to_position(index_in_vec, nrow(dat), ncol(dat))
  index_out_mat <- .convert_index_to_position(index_out_vec, nrow(dat), ncol(dat))
  n_in <- length(index_in_vec); n_out <- length(index_out_vec)

  if(row) z = 1 else z = 2
  func <- function(x){
    index_in <- index_in_mat[which(index_in_mat[,z] == x), -z+3]
    index_out <- index_out_mat[which(index_out_mat[,z] == x), -z+3]

    .estimate_row(dat, initial_mat[x,], latent_mat, x, index_in, index_out, n_in, n_out,
                  tol, max_iter, row, verbose = F)
  }

  res <- foreach::"%dopar%"(foreach::foreach(x = 1:nrow(initial_mat)), func(x))

  do.call(rbind, res)
}

.convert_index_to_position <- function(vec, num_row, num_col){
  stopifnot(all(vec <= num_col * num_row), all(vec %% 1 == 0), all(vec >= 1),
            length(unique(vec)) == length(vec))

  tmp <- vec %% num_row
  tmp[tmp == 0] <- num_row
  cbind(tmp, ceiling(vec / num_row))
}

# verbose only for debugging purposes only
.estimate_row <- function(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out,
                          n_in, n_out, tol = 1e-5, max_iter = 500, row = T, verbose = F){
  vec <- initial_vec
  obj_old <- .evaluate_objective_single(dat, vec, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row)

  if(verbose) {obj_vec <- rep(NA, max_iter+1); obj_vec[1] <- obj_old}
  subgrad <- .subgradient_vec(dat, vec, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row)
  k <- .initialize_k(dat, vec, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row, subgrad, obj_old)

  while(TRUE){
    subgrad <- .subgradient_vec(dat, vec, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row)
    vec <- vec - 1/k*subgrad
    obj_new <- .evaluate_objective_single(dat, vec, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row)

    if(abs(obj_old - obj_new) < tol) break()
    if(k > max_iter) break()

    if(verbose) obj_vec[k+1] <- obj_new

    obj_old <- obj_new
    k <- k+1
  }

  if(verbose) list(vec = vec, obj_vec = obj_vec) else vec
}

.subgradient_vec <- function(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row = T){
  stopifnot(ncol(latent_mat) == length(initial_vec), (ncol(dat) == nrow(latent_mat) |
                                                        nrow(dat) == nrow(latent_mat)))
  stopifnot((row & fixed_idx <= nrow(dat)) | (!row & fixed_idx <= ncol(dat)))
  if((length(index_in) > 0 | length(index_out) > 0)) stopifnot(max(c(index_in, index_out)) <= prod(dim(dat)))
  stopifnot(length(intersect(index_in, index_out)) == 0)
  stopifnot(length(index_in) <= n_in, length(index_out) <= n_out)

  if(length(index_out) > 0){
    vec <- initial_vec %*% t(latent_mat[index_out,,drop = F])
    vec <- sapply(vec, function(x){max(0, x)})
    second_term <- vec%*% latent_mat[index_out,,drop = F]
  } else second_term <- rep(0, length(initial_vec))

  if(length(index_in) > 0){
    if(row) tmp <- dat[fixed_idx,index_in] else tmp <- dat[index_in,fixed_idx]
    first_term <- -(tmp - initial_vec %*% t(latent_mat[index_in,,drop = F])) %*% latent_mat[index_in,,drop = F]
  } else first_term <- rep(0, length(initial_vec))

  as.numeric(first_term/n_in + second_term/n_out)
}

.evaluate_objective_single <- function(dat, initial_vec, latent_mat, fixed_idx, index_in,
                                       index_out, n_in, n_out, row = T){
  stopifnot(ncol(latent_mat) == length(initial_vec))
  stopifnot((row & ncol(dat) == nrow(latent_mat)) | (!row & nrow(dat) == nrow(latent_mat)))
  stopifnot(length(fixed_idx) <= nrow(dat))
  if((length(index_in) > 0 | length(index_out) > 0)) stopifnot(max(c(index_in, index_out)) <= prod(dim(dat)))
  stopifnot(length(intersect(index_in, index_out)) == 0)
  stopifnot(length(index_in) <= n_in, length(index_out) <= n_out)

  if(length(index_out) > 0){
    vec <- initial_vec %*% t(latent_mat[index_out,])
    vec <- sapply(vec, function(x){max(0, x)})
    second_term <- sum(vec^2)
  } else second_term <- 0

  if(length(index_in) > 0){
    if(row) tmp <- dat[fixed_idx,index_in] else tmp <- dat[index_in,fixed_idx]
    first_term <- sum((tmp - initial_vec %*% t(latent_mat[index_in,,drop = F]))^2)
  } else first_term <- 0

  as.numeric(first_term/(2*n_in) + second_term/(2*n_out))
}

.evaluate_objective_full <- function(dat, u_mat, v_mat, index_in, index_out){
  pred_mat <- u_mat %*% t(v_mat)
  as.numeric(sum((dat[index_in] - pred_mat[index_in])^2)/length(index_in) +
    sum(pmax(0, pred_mat[index_out])^2)/length(index_out))
}

.initialize_k <- function(dat, vec, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row, subgrad, obj_old){
  k <- 1

  seq_vec <- seq(1, 1000, by = 10)

  while(TRUE){
    vec2 <- vec - 1/seq_vec[k]*subgrad
    obj_new <- .evaluate_objective_single(dat, vec2, latent_mat, fixed_idx, index_in, index_out, n_in, n_out, row)
    if(obj_new < obj_old) break()

    if(k >= length(seq_vec)) break()
    k <- k+1
  }

  seq_vec[k]
}
