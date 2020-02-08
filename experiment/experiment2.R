rm(list=ls())
set.seed(1)
n <- 50
u_mat <- matrix(abs(rnorm(50)), ncol = 1)
v_mat <- matrix(abs(rnorm(50)), ncol = 1)
pred_mat <- u_mat %*% t(v_mat)

dat <- pred_mat

for(i in 1:10){
  for(j in 1:4){
    dat[i,j] <- abs(stats::rnorm(1, mean = 1/pred_mat[i,j], sd = 1/(2*pred_mat[i,j])))
  }
}

family = "curved_gaussian"
k = 2
tol = 1e-3
max_val = 100
max_iter = 10
verbose = F

direction <- .dictate_direction(family)
if(!is.na(direction) & !is.na(max_val)){
  stopifnot((direction == ">=" & max_val > 0) | (direction == "<=" & max_val < 0))
}

# initialize
dat <- .matrix_completion(dat, k = k)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

min_val <- min(dat[which(dat > 0)])
dat[which(dat <= 0)] <- min_val/2
pred_mat <- .mean_transformation(dat, family)
direction <- .dictate_direction(family)

n <- nrow(dat); d <- ncol(dat)
# pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)

#####################

min_val <- max(min(dat[which(dat > 0)]), tol)
dat[which(dat <= min_val)] <- min_val/2
pred_mat <- .mean_transformation(dat, family, ...)
direction <- .dictate_direction(family)

# .project_rank_feasibility(pred_mat, k = k, direction = direction,
#                           max_val = max_val)$matrix

####################

mat <- pred_mat
stopifnot(!is.na(max_val) | !is.na(direction))
if(!is.na(max_val) & !is.na(direction)) stopifnot((direction == "<=" & max_val < 0) | (direction == ">=" & max_val > 0))

mat_org <- mat
iter <- 1

while(iter < max_iter){
  res <- .svd_projection(mat, k = k, factors = T)
  mat <- res$u_mat %*% t(res$v_mat)

  if(is.na(direction)){
    if(all(abs(mat) <= max_val+tol)) return(list(matrix = mat, iter = iter))
  } else if (direction == "<=") {
    if(all(mat < 0) && (is.na(max_val) || all(mat > max_val-tol))) return(list(matrix = mat, iter = iter))
  } else {
    if(all(mat > 0) && (is.na(max_val) || all(mat < max_val+tol))) return(list(matrix = mat, iter = iter))
  }

  mat <- .absolute_threshold(mat, direction, max_val)

  iter <- iter + 1
}

mat <- .absolute_threshold(mat_org, direction, max_val)
# res <- .dcsbm_projection(mat, k)

#########################

res <- .spectral_clustering(mat, k)
row_clustering <- res$row_clustering
col_clustering <- res$col_clustering

res <- .form_prediction_dcsbm(mat, row_clustering, col_clustering)

base_mat <- sapply(1:ncol(mat), function(i){
  res$o_mat[row_clustering, col_clustering[i]]
})

base_mat <- diag(res$theta_row) %*% base_mat %*% diag(res$theta_col)
