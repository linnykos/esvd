set.seed(10)
dat <- matrix(rnbinom(40, size = 10, prob = 0.5), nrow = 5, ncol = 5)
cv_trials <- 3
missing_idx_list <- lapply(1:cv_trials, function(i){
  construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 1)
})

init <- initialization(dat, family = "neg_binom", max_val = 100, scalar = 10)
fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                  max_iter = 10, max_val = 100,
                  family = "neg_binom", scalar = 10)

#############################
u_mat = init$u_mat
v_mat = init$v_mat
max_iter = 10
max_val = 100
family = "neg_binom"
scalar = 10
verbose = T
return_path = F
ncores = NA

if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)
stopifnot(length(which(dat > 0)) > 0)
stopifnot(length(which(dat < 0)) == 0)
stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
          ncol(u_mat) == ncol(v_mat))
attr(dat, "family") <- family
if(!is.na(max_val)) stopifnot(max_val >= 0)

tmp <- .check_rank(u_mat, v_mat)
u_mat <- tmp$u_mat; v_mat <- tmp$v_mat
k <- ncol(u_mat); n <- nrow(u_mat); p <- nrow(v_mat)

idx <- which(!is.na(dat))
min_val <- min(dat[which(dat > 0)])
dat[which(dat == 0)] <- min_val/2
direction <- .dictate_direction(attr(dat, "family"))
if(!is.na(direction) && direction == "<=" && !is.na(max_val)) max_val <- -max_val

current_obj <- Inf
next_obj <- .evaluate_objective(dat, u_mat, v_mat, scalar = scalar)
obj_vec <- c(next_obj)
if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))
if(return_path) res_list <- list(list(u_mat = u_mat, v_mat = v_mat)) else res_list <- NA

#############

current_obj <- next_obj

nat_mat <- u_mat%*%t(v_mat)
if(!is.na(direction)){ if(direction == "<=") { stopifnot(all(nat_mat <= 0))
} else { stopifnot(all(nat_mat >= 0)) }
}

# reparameterize
tmp <- .reparameterize(u_mat, v_mat, equal_covariance = F)
u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

# u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, parallelized = !is.na(ncores), scalar = scalar)

##############
current_mat = u_mat
other_mat = v_mat
left = T
max_val = max_val

stopifnot(ncol(current_mat) == ncol(other_mat))
if(left) {
  stopifnot(ncol(dat) == nrow(other_mat))
} else {
  stopifnot(nrow(dat) == nrow(other_mat))
}


i = 1
if(left) {
  dat_vec <- dat[i,]
} else {
  dat_vec <- dat[,i]
}
attr(dat_vec, "family") <- attr(dat, "family")
# if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
#                                                           other_mat, max_val = max_val,
#                                                           n = n, p = p, scalar = scalar)

###########
current_vec = current_mat[i,]
stopifnot(length(which(!is.na(dat_vec))) > 0, length(dat_vec) == nrow(other_mat))

direction <- .dictate_direction(attr(dat_vec, "family"))
current_obj <- Inf
next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat, n = n, p = p, scalar = scalar)
iter <- 1
current_obj <- next_obj

grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat, n = n, p = p, scalar = scalar)
# step_vec <- .frank_wolfe(grad_vec, other_mat,
#                          direction = direction, other_bound = max_val)

################
direction = direction
other_bound = max_val
tol = 0.001

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

res <- suppressWarnings(clplite::clp_solve(objective_in, constr_mat, constr_lb, constr_ub, var_lb, var_ub, max = F))


