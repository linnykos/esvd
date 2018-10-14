rm(list=ls())
source("../experiment/Week25_simulation_generator.R")
library(singlecell)

set.seed(10)
simulation <- .data_generator_exponential(total = 200, col_drop = F)
dat <- simulation$dat

zz <- dat[dat > 0]
max_val <- -1/mean(zz[zz < quantile(zz, probs = 0.2)])
max_iter <- 50

set.seed(10)
dat <- simulation$obs_mat
u_mat <- simulation$cell_mat
v_mat <- simulation$gene_mat
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
next_obj <- singlecell:::.evaluate_objective.exponential(dat, u_mat, v_mat)
obj_vec <- c(next_obj)
if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))

#####################
# do one iteration
current_obj <- next_obj

u_mat <- singlecell:::.optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, !is.na(cores))
print("U done")
v_mat <- singlecell:::.optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, !is.na(cores))

next_obj <- singlecell:::.evaluate_objective.exponential(dat, u_mat, v_mat)

if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))

obj_vec <- c(obj_vec, next_obj)

######################
# zoom in on the next iteration
current_mat <- u_mat
other_mat <- v_mat
left <- T
paallelized <- F

# do as many iterations that don't break
for(i in 1:138){
  if(left) { dat_vec <- dat[i,] } else { dat_vec <- dat[,i] }
  class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
  if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                            other_mat, max_val = max_val)
}

#######################
# focus on the one iteration that breaks
i <- 139
if(left) { dat_vec <- dat[i,] } else { dat_vec <- dat[,i] }
class(dat_vec) <- c(class(dat)[1], class(dat_vec)[length(class(dat_vec))])
if(any(!is.na(dat_vec))) current_mat[i,] <- .optimize_row(dat_vec, current_mat[i,],
                                                          other_mat, max_val = max_val)

# check to see if in bounds before starting
min(other_mat %*% current_mat[i,]) >= max_val #yes

current_vec <- current_mat[i,]
max_iter <- 100

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

while(abs(current_obj - next_obj) > 1e-6 & iter < 8){
  current_obj <- next_obj

  grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat)
  step_vec <- .frank_wolfe(grad_vec, other_mat, which(!is.na(dat_vec)),
                           direction = direction, other_bound = max_val)
  step_size <- .binary_search(dat_vec, current_vec, step_vec, other_mat)
  current_vec <- (1-step_size)*current_vec + step_size*step_vec

  next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)

  iter <- iter + 1
}

######################
# focus on the iteration that breaks: iter = 8
current_obj <- next_obj

#check current status
min(other_mat %*% current_vec) >= max_val #yes

grad_vec <- .gradient_vec(dat_vec, current_vec, other_mat)
# step_vec <- .frank_wolfe(grad_vec, other_mat, which(!is.na(dat_vec)),
#                          direction = direction, other_bound = max_val)
# step_size <- .binary_search(dat_vec, current_vec, step_vec, other_mat)
# current_vec <- (1-step_size)*current_vec + step_size*step_vec
# next_obj <- .evaluate_objective_single(dat_vec, current_vec, other_mat)

#######################
# start diving into frank_wolfe
idx <- which(!is.na(dat_vec))
tol = 1e-6
direction = "<="
other_bound = max_val

k <- length(grad_vec)
other_mat <- other_mat[idx,,drop = F]
other_direction <- ifelse(direction == "<=", ">=", "<=")

objective_in <- c(grad_vec, -grad_vec)

const_mat <- t(apply(other_mat, 1, function(v){ c(v, -v) }))
if(!is.na(other_bound)) const_mat <- rbind(const_mat, const_mat)

const_dir <- rep(direction, nrow(other_mat))
if(!is.na(other_bound)) const_dir <- c(const_dir, rep(other_direction, nrow(other_mat)))

if(direction == "<=") s <- -1 else s <- 1
const_rhs <- rep(s*tol, nrow(other_mat))
if(!is.na(other_bound)) const_rhs <- c(const_rhs, rep(other_bound, nrow(other_mat)))

res <- lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)
# Error: no feasible solution found ???

########################
# investigation
range(const_mat %*% c(current_vec, 0, 0))

zz = const_mat %*% c(current_vec, tol, tol)
all(zz[1:nrow(other_mat)] == zz[(nrow(other_mat)+1):(2*nrow(other_mat))])
all(zz[1:nrow(other_mat)] <= const_rhs[1:nrow(other_mat)]) &
  all(zz[(nrow(other_mat)+1):(2*nrow(other_mat))] >= const_rhs[(nrow(other_mat)+1):(2*nrow(other_mat))])
all(const_dir[1:nrow(other_mat)] == "<=")
all(const_dir[(nrow(other_mat)+1):(2*nrow(other_mat))] == ">=")
all(c(current_vec, tol, tol) >= 0)

## try removing random rows to see what's going on
val <- 240 #val of 182 would ruin it...
idx <- c(1:val, 241:(240+val))
const_rhs2 = const_rhs;const_rhs2[241:480] <- -max_val
const_dir2 <- rep("<=", 480)
const_mat2 <- const_mat; const_mat2[241:480,] <- -const_mat2[241:480,]
res <- lpSolve::lp("min", objective_in, const_mat2[idx,],
                   const_dir2[idx], const_rhs2[idx])
res_yy <- res$solution[1:2] - res$solution[3:4]
yy <- const_mat %*% res$solution
range(yy)
all(yy[1:nrow(other_mat)] <= const_rhs[1:nrow(other_mat)]) &
all(yy[(nrow(other_mat)+1):(2*nrow(other_mat))] >= const_rhs[(nrow(other_mat)+1):(2*nrow(other_mat))])

# try a different strategy
for(i in 1:181){
  res <- lpSolve::lp("min", objective_in, const_mat[c(i,182,240+i,240+182),],
                     const_dir[c(i,182,240+i,240+182)], const_rhs[c(i,182,240+i,240+182)])
  if(res$status != 0) stop()
}

## try a different optimizaiton program
objective_in <- c(grad_vec, -grad_vec, rep(0, length(grad_vec)))

const_mat <- t(apply(other_mat, 1, function(v){ c(v, -v, rep(0, length(v))) }))
const_mat <- rbind(const_mat, cbind(diag(2), -diag(2), -diag(2))) #x_+ - x_- <= s
const_mat <- rbind(const_mat, cbind(-diag(2), diag(2), -diag(2))) #-x_+ + x_- <= s
const_mat <- rbind(const_mat, cbind(matrix(0, ncol = 4, nrow = 2), diag(2)))

const_dir <- rep("<=", nrow(other_mat))
const_dir <- c(const_dir, rep("<=", 2), rep("<=", 2), rep("<=", 2))

const_rhs <- rep(-tol, nrow(other_mat))
const_rhs <- c(const_rhs, rep(0, 4), rep(15, 2))

lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)

##another attempt
objective_in <- c(grad_vec, -grad_vec)

const_mat <- t(apply(other_mat, 1, function(v){ c(v, -v) }))
const_mat <- rbind(const_mat, cbind(diag(2), matrix(0,2,2)))
const_mat <- rbind(const_mat, cbind(matrix(0,2,2), diag(2)))

const_dir <- rep("<=", nrow(other_mat))
const_dir <- c(const_dir, rep("<=", 2), rep("<=", 2))

const_rhs <- rep(-tol, nrow(other_mat))
const_rhs <- c(const_rhs, rep(15, 4))

lpSolve::lp("min", objective_in, const_mat, const_dir, const_rhs)

#####################

# let's try using clplite
objective_in <- c(grad_vec)

const_mat <- other_mat
const_lb <- rep(max_val, nrow(const_mat))
const_ub <- rep(-tol, nrow(const_mat))
var_lb <- rep(-abs(max_val), nrow(const_mat))
var_ub <- rep(abs(max_val), nrow(const_mat))

res <- clplite::clp_solve(objective_in, const_mat, const_lb, const_ub, var_lb, var_ub, max = F)
