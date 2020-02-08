rm(list=ls())
set.seed(10)
dat <- matrix(rexp(100), nrow = 10, ncol = 10)

dat2 <- dat; dat2[2,1] <- NA
dat3 <- dat
for(i in 1:nrow(dat3)){
  dat3[i, sample(1:ncol(dat3), 2)] <- NA
}

init <- initialization(dat, max_val = 100, family = "curved_gaussian")

range(init$u_mat %*% t(init$v_mat))

# res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
#                          max_val = 100, max_iter = 100, family = "curved_gaussian")

##################
u_mat = init$u_mat
v_mat = init$v_mat
max_val = 100
max_iter = 100
family = "curved_gaussian"

reparameterize = T
tol = 1e-3
verbose = F
return_path = F
cores = NA

if(!is.na(cores)) doMC::registerDoMC(cores = cores)
stopifnot(length(which(dat > 0)) > 0)
stopifnot(length(which(dat < 0)) == 0)
stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
          ncol(u_mat) == ncol(v_mat))
k <- ncol(u_mat)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

idx <- which(!is.na(dat))
min_val <- min(dat[which(dat > 0)])
dat[which(dat == 0)] <- min_val/2
direction <- .dictate_direction(class(dat)[1])

current_obj <- Inf
next_obj <- .evaluate_objective(dat, u_mat, v_mat)
obj_vec <- c(next_obj)
if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))
if(return_path) res_list <- list(list(u_mat = u_mat, v_mat = v_mat)) else res_list <- NA


#while((is.na(tol) | abs(current_obj - next_obj) > tol) & length(obj_vec) < max_iter){
  print(length(obj_vec))
  print(range(u_mat %*% t(v_mat)))
  current_obj <- next_obj

  pred_mat <- u_mat%*%t(v_mat)
  if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(pred_mat <= 0))
  } else { stopifnot(all(pred_mat >= 0)) }
  }

  if(reparameterize){
    svd_res <- .svd_projection(pred_mat, k = k, factors = T)
    u_mat <- svd_res$u_mat; v_mat <- svd_res$v_mat
  }

  print(range(u_mat %*% t(v_mat)))

  u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, parallelized = !is.na(cores))

  print(range(u_mat %*% t(v_mat)))

  pred_mat <- u_mat%*%t(v_mat)
  if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(pred_mat <= 0))
  } else { stopifnot(all(pred_mat >= 0)) }
  }

  print(range(u_mat %*% t(v_mat)))

  if(reparameterize){
    svd_res <- .svd_projection(pred_mat, k = k, factors = T)
    u_mat <- svd_res$u_mat; v_mat <- svd_res$v_mat
  }

  print(paste0("? ", range(u_mat %*% t(v_mat))))

  v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, parallelized = !is.na(cores))

  print(range(u_mat %*% t(v_mat)))

  next_obj <- .evaluate_objective(dat, u_mat, v_mat)


  if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))
  if(return_path) res_list[[length(res_list)+1]] <- list(u_mat = u_mat, v_mat = v_mat)

  obj_vec <- c(obj_vec, next_obj)
#}

  print(length(obj_vec))
  print(range(u_mat %*% t(v_mat)))
  current_obj <- next_obj

  pred_mat <- u_mat%*%t(v_mat)
  if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(pred_mat <= 0))
  } else { stopifnot(all(pred_mat >= 0)) }
  }

  if(reparameterize){
    svd_res <- .svd_projection(pred_mat, k = k, factors = T)
    u_mat <- svd_res$u_mat; v_mat <- svd_res$v_mat
  }

  print(range(u_mat %*% t(v_mat)))

  u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, parallelized = !is.na(cores))

  print(range(u_mat %*% t(v_mat)))

  pred_mat <- u_mat%*%t(v_mat)
  if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(pred_mat <= 0))
  } else { stopifnot(all(pred_mat >= 0)) }
  }

  print(range(u_mat %*% t(v_mat)))


mat = u_mat %*% t(v_mat)
svd_res <- .svd_projection(mat, k = k, factors = T)
print(paste0("? ", range(svd_res$u_mat %*% t(svd_res$v_mat))))
