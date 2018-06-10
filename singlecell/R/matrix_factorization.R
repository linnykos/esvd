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

  if(row) z = 1 else z = 2
  func <- function(x){
    index_in <- index_in_mat[which(index_in_mat[,z] == x), -z+3]
    index_out <- index_out_mat[which(index_out_mat[,z] == x), -z+3]

    .estimate_row(dat, initial_mat[x,], latent_mat, x, index_in, index_out,
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
                          tol = 1e-5, max_iter = 500, row = T, verbose = F){
  k <- 1
  vec <- initial_vec
  obj_old <- .evaluate_objective_single(dat, vec, latent_mat, fixed_idx, index_in, index_out, row)

  if(verbose) {obj_vec <- rep(NA, max_iter+1); obj_vec[1] <- obj_old}

  while(TRUE){
    subgrad <- .subgradient_vec(dat, vec, latent_mat, fixed_idx, index_in, index_out, row)
    vec <- vec - 1/k*subgrad
    obj_new <- .evaluate_objective_single(dat, vec, latent_mat, fixed_idx, index_in, index_out, row)

    if(abs(obj_old - obj_new) < tol) break()
    if(k > max_iter) break()

    if(verbose) obj_vec[k+1] <- obj_new

    obj_old <- obj_new
    k <- k+1
  }

  if(verbose) list(vec = vec, obj_vec = obj_vec) else vec
}

.subgradient_vec <- function(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out, row = T){
  stopifnot(ncol(latent_mat) == length(initial_vec), (ncol(dat) == nrow(latent_mat) |
                                                        nrow(dat) == nrow(latent_mat)))
  stopifnot((row & fixed_idx <= nrow(dat)) | (!row & fixed_idx <= ncol(dat)))
  stopifnot(max(c(index_in, index_out)) <= prod(dim(dat)))
  stopifnot(length(intersect(index_in, index_out)) == 0)

  vec <- initial_vec %*% t(latent_mat[index_out,])
  vec <- sapply(vec, function(x){max(0, x)})
  second_term <- 2*vec%*% latent_mat[index_out,]

  if(row) tmp <- dat[fixed_idx,index_in] else tmp <- dat[index_in,fixed_idx]
  -(tmp - initial_vec %*% t(latent_mat[index_in,])) %*% latent_mat[index_in,]/length(index_in) +
    second_term/length(index_out)
}

.evaluate_objective_single <- function(dat, initial_vec, latent_mat, fixed_idx, index_in,
                                       index_out, row = T){
  stopifnot(ncol(latent_mat) == length(initial_vec))
  stopifnot((row & ncol(dat) == nrow(latent_mat)) | (!row & nrow(dat) == nrow(latent_mat)))
  stopifnot(length(fixed_idx) <= nrow(dat))
  stopifnot(max(c(index_in, index_out)) <= prod(dim(dat)))
  stopifnot(length(intersect(index_in, index_out)) == 0)

  vec <- initial_vec %*% t(latent_mat[index_out,])
  vec <- sapply(vec, function(x){max(0, x)})
  second_term <- sum(vec^2)

  if(row) tmp <- dat[fixed_idx,index_in] else tmp <- dat[index_in,fixed_idx]
  sum((tmp - initial_vec %*% t(latent_mat[index_in,]))^2)/(2*length(index_in)) +
    second_term/(2*length(index_out))
}

.evaluate_objective_full <- function(dat, u_mat, v_mat, index_in, index_out){
  pred_mat <- u_mat %*% t(v_mat)
  sum((dat[index_in] - pred_mat[index_in])^2)/length(index_in) +
    sum(pmax(0, pred_mat[index_out])^2)/length(index_out)
}
