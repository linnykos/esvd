rm(list=ls())
set.seed(10)
u_mat <- cbind(c(rep(0.1, 5), rep(0.5, 5)), c(rep(0.3, 5), rep(1, 5)))
v_mat <- cbind(c(rep(0.5, 2), rep(0.2, 2)), c(rep(0.1, 2), rep(1, 2)))
s_vec <- c(5:15)

dat <- matrix(0, nrow = 10, ncol = 4)
for(i in 1:10){
  for(j in 1:4){
    dat[i,j] <- stats::rpois(1, lambda = s_vec[i]*u_mat[i,]%*%v_mat[j,])
  }
}

extra_weights <- rowSums(dat)
init <- .initialization(dat, family = "poisson", max_val = 100,
                        extra_weights = extra_weights)

###########

u_mat = init$u_mat
v_mat = init$v_mat
max_iter = 5
max_val = 100
family = "poisson"
reparameterize = F
cores = NA
verbose = F
tol = 1e-3

if(!is.na(cores)) doMC::registerDoMC(cores = cores)
stopifnot(length(which(dat > 0)) > 0)
stopifnot(length(which(dat < 0)) == 0)
stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
          ncol(u_mat) == ncol(v_mat))
k <- ncol(u_mat)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

idx <- which(dat == 0)
min_val <- min(dat[which(dat > 0)])
dat[which(dat == 0)] <- min_val/2

current_obj <- Inf
next_obj <- .evaluate_objective(dat, u_mat, v_mat)
obj_vec <- c(next_obj)
if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))
current_obj <- next_obj

# u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, extra_weights = extra_weights,
#                        !is.na(cores))
##############
current_mat = u_mat
other_mat = v_mat
left = T
parallelized = F

stopifnot(length(class(dat)) == 2)

stopifnot(ncol(current_mat) == ncol(other_mat))
if(left) {
  stopifnot(ncol(dat) == nrow(other_mat))
} else {
  stopifnot(nrow(dat) == nrow(other_mat))
}

for(i in 1:nrow(dat)){
  print(i)
  if(left) {
    dat_vec <- dat[i,]; other_mat <- other_mat * extra_weights[i]
  } else {
    dat_vec <- dat[,i]; other_mat <- diag(extra_weights) %*% other_mat
  }
  class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
  if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                            other_mat, max_val = max_val)
}
