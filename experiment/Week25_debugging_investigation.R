rm(list=ls())
load("../experiment/Week25_simulation_exponential2.RData")

#investigate dat_impute
# res_withimpute <- singlecell:::.fit_factorization(dat_impute, init_impute$u_mat, init_impute$v_mat,
#                                                   verbose = T, family = "exponential",
#                                                   max_iter = max_iter,
#                                                   cores = 15, max_val = max_val)

u_mat <- init_impute$u_mat
v_mat <- init_impute$v_mat
verbose <- T
family <- "exponential"
cores <- NA
tol <- 1e-3

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green

plot(u_mat[,1], u_mat[,2],
     xlim = range(c(u_mat[,1], 0)),
     ylim = range(c(u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")

############

dat <- dat_impute
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])
idx <- which(dat == 0)
min_val <- min(dat[which(dat > 0)])
dat[which(dat == 0)] <- min_val/2

current_obj <- Inf
next_obj <- .evaluate_objective(dat, u_mat, v_mat)
obj_vec <- c(next_obj)

while(current_obj - next_obj > tol & length(obj_vec) < max_iter){
  current_obj <- next_obj

  u_mat2 <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, !is.na(cores))
  obj_tmp1 <- .evaluate_objective(dat, u_mat2, v_mat)

  if(obj_tmp1 > current_obj) stop()

  v_mat2 <- .optimize_mat(dat, v_mat, u_mat2, left = F, max_val = max_val, !is.na(cores))
  next_obj <- .evaluate_objective(dat, u_mat2, v_mat2)

  if(next_obj > obj_tmp1) stop()

  u_mat <- u_mat2
  v_mat <- v_mat2

  if(verbose) print(paste0("Iter ", length(obj_vec), ": Decrease is ", current_obj - next_obj))

  obj_vec <- c(obj_vec, next_obj)
}
