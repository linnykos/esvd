rm(list=ls())
load("../results/step2_naive_svd_spca.RData")

max_val <- 5000
scalar_vec <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100)
res_list <- vector("list", length(scalar_vec))

i <- 10
print(paste0("Trying scalar value = ", scalar_vec[i]))
init <- singlecell::initialization(dat_impute_NA, family = family,
                                   k = k, max_val = max_val)

stopifnot(all(init$u_mat %*% t(init$v_mat) > 0))
range(init$u_mat %*% t(init$v_mat) )

res_list[[i]] <- singlecell::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                               family = family, reparameterize = T,
                                               max_iter = 25, max_val = NA,
                                               scalar = scalar_vec[i],
                                               return_path = F, cores = NA,
                                               verbose = T)
##########################

dat = dat_impute_NA
u_mat = init$u_mat
v_mat = init$v_mat
max_iter = 25
reparameterize = T
max_val = NA
scalar = scalar_vec[i]
return_path = F
cores = NA
verbose = T
tol = 1e-3

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
next_obj <- .evaluate_objective(dat, u_mat, v_mat, scalar = scalar)
obj_vec <- c(next_obj)
if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))
if(return_path) res_list <- list(list(u_mat = u_mat, v_mat = v_mat)) else res_list <- NA

current_obj <- next_obj

pred_mat <- u_mat%*%t(v_mat)
if(direction == "<=") stopifnot(all(pred_mat < 0)) else  stopifnot(all(pred_mat >= 0))

if(reparameterize){
  svd_res <- .svd_projection(pred_mat, k = k, factors = T)
  u_mat <- svd_res$u_mat; v_mat <- svd_res$v_mat
}

#######################

left = T
parallelized = !is.na(cores)
current_mat <- u_mat
other_mat <- v_mat

stopifnot(length(class(dat)) == 2)

stopifnot(ncol(current_mat) == ncol(other_mat))
if(left) {
  stopifnot(ncol(dat) == nrow(other_mat))
} else {
  stopifnot(nrow(dat) == nrow(other_mat))
}

i = 93


if(left) {
  dat_vec <- dat[i,]
} else {
  dat_vec <- dat[,i]
}

current_vec <- current_mat[i,]

class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                          other_mat, max_val = max_val,
                                                          scalar = scalar)

range(current_mat[i,]%*%t(other_mat))

#################

max_iter = 100

stopifnot(length(which(!is.na(dat_vec))) > 0)

direction <- .dictate_direction(class(dat_vec)[1])
current_obj <- Inf
next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat,
                                       scalar = scalar)
iter <- 1
while(abs(current_obj - next_obj) > 1e-6 & iter < 2){
  current_obj <- next_obj

  grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat, scalar = scalar)
  step_vec <- .frank_wolfe(grad_vec, other_mat,
                           direction = direction, other_bound = max_val)
  step_size <- .binary_search(dat_vec, current_vec, step_vec, other_mat, scalar = scalar)
  current_vec <- (1-step_size)*current_vec + step_size*step_vec

  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat, scalar = scalar)

  if(!is.na(max_val)) stopifnot(all(abs(other_mat %*% current_vec) <= abs(max_val)+1e-6))
  iter <- iter + 1
}

current_obj <- next_obj

grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat, scalar = scalar)
range(current_vec%*%t(other_mat))

# step_vec <- .frank_wolfe(grad_vec, other_mat,
#                          direction = direction, other_bound = max_val)

#############

tol = 0.001
other_bound = NA

k <- length(grad_vec)
other_direction <- ifelse(direction == "<=", ">=", "<=")
objective_in <- grad_vec
constr_mat <- other_mat
k <- nrow(constr_mat)

if(direction == "<="){
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

###############

zz <- res$solution
zz <- current_vec

all(constr_mat %*% zz >= constr_lb)
# min(constr_mat %*% zz)
all(constr_mat %*% zz <= constr_ub)
all(zz >= var_lb)
all(zz <= var_ub)


