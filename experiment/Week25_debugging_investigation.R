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

while((is.na(tol) | abs(current_obj - next_obj) > tol) & length(obj_vec) < max_iter){
  current_obj <- next_obj

  u_mat <- singlecell:::.optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, !is.na(cores))
  print("U done")
  v_mat <- singlecell:::.optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, !is.na(cores))

  next_obj <- singlecell:::.evaluate_objective.exponential(dat, u_mat, v_mat)

  if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))

  obj_vec <- c(obj_vec, next_obj)
}
