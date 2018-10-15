rm(list=ls())
source("../experiment/Week25_simulation_generator.R")
load("../experiment/Week25_simulation_exponential.RData")

set.seed(10)
init_impute <- singlecell:::.initialization(dat_impute, family = "exponential",
                                            max_val = max_val)

###########

dat <- dat2
u_mat <- init_impute$u_mat
v_mat <- init_impute$v_mat
verbose = T
family = "exponential"
tol = NA
cores = NA

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
u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, !is.na(cores))
# v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, !is.na(cores))

#####################

current_mat = v_mat
other_mat = u_mat
left = F
parallelized = !is.na(cores)

for(i in 1:14){
  if(left) { dat_vec <- dat[i,] } else { dat_vec <- dat[,i] }
  class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
  if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                            other_mat, max_val = max_val)
}

i = 15
if(left) { dat_vec <- dat[i,] } else { dat_vec <- dat[,i] }
class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])

#################
range(other_mat %*% current_mat[i,])

max_iter = 100
current_vec = current_mat[i,]
if(class(dat_vec)[1] == "exponential"){
  direction = "<="
} else if(class(dat_vec)[1] == "gaussian"){
  direction = ">="
} else {
  stop("input vector in .optimize_row() does not have a proper class")
}

current_obj <- Inf
next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)
iter <- 1

while(abs(current_obj - next_obj) > 1e-6 & iter < 2){
  current_obj <- next_obj

  grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat)
  step_vec <- .frank_wolfe(grad_vec, other_mat, which(!is.na(dat_vec)),
                           direction = direction, other_bound = max_val)
  step_size <- .binary_search(dat_vec, current_vec, step_vec, other_mat)
  current_vec <- (1-step_size)*current_vec + step_size*step_vec

  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)

  iter <- iter + 1
}

grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat)
step_vec <- .frank_wolfe(grad_vec, other_mat, which(!is.na(dat_vec)),
                         direction = direction, other_bound = max_val)

range(other_mat[which(!is.na(dat_vec)),] %*% step_vec) # BUG
#############
