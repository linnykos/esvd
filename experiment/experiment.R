rm(list=ls())
set.seed(10)
dat <- matrix(rnbinom(40, size = 10, prob = 0.5), nrow = 5, ncol = 5)
class(dat) <- c("neg_binom", class(dat)[length(class(dat))])
init <- initialization(dat, family = "neg_binom", max_val = -100, r = 10)

u_mat = init$u_mat
v_mat = init$v_mat
max_iter = 5
max_val = -100
family = "neg_binom"
r = 10
reparameterize = T
tol = 1e-3
verbose = F
return_path = F
cores = NA

###############

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
next_obj <- .evaluate_objective(dat, u_mat, v_mat, r=r)
obj_vec <- c(next_obj)
if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))
if(return_path) res_list <- list(list(u_mat = u_mat, v_mat = v_mat)) else res_list <- NA


while((is.na(tol) | abs(current_obj - next_obj) > tol) & length(obj_vec) < max_iter){
  current_obj <- next_obj

  pred_mat <- u_mat%*%t(v_mat)
  if(direction == "<=") {
    stopifnot(all(pred_mat <= 0))
  } else if(direction == ">=") {
    stopifnot(all(pred_mat >= 0))
  }

  if(reparameterize){
    svd_res <- .svd_projection(pred_mat, k = k, factors = T)
    u_mat <- svd_res$u_mat; v_mat <- svd_res$v_mat
  }

  u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, parallelized = !is.na(cores), r=r)

  pred_mat <- u_mat%*%t(v_mat)
  if(direction == "<=") {
    stopifnot(all(pred_mat <= 0))
  } else if(direction == ">=") {
    stopifnot(all(pred_mat >= 0))
  }

  if(reparameterize){
    svd_res <- .svd_projection(pred_mat, k = k, factors = T)
    u_mat <- svd_res$u_mat; v_mat <- svd_res$v_mat
  }

  v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, parallelized = !is.na(cores), r=r)

  next_obj <- .evaluate_objective(dat, u_mat, v_mat, r=r)


  if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))
  if(return_path) res_list[[length(res_list)+1]] <- list(u_mat = u_mat, v_mat = v_mat)

  obj_vec <- c(obj_vec, next_obj)
}

tmp <- .reparameterize(u_mat, v_mat)
u_mat <- tmp$u_mat; v_mat <- tmp$v_mat
pred_mat <- u_mat%*%t(v_mat)
if(direction == "<=") {
  stopifnot(all(pred_mat < 0))
} else if(direction == ">=") {
  stopifnot(all(pred_mat >= 0))
}

list(u_mat = u_mat, v_mat = v_mat, obj_vec = obj_vec, res_list = res_list)
