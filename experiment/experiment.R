rm(list=ls())

x = 3
set.seed(13*x)
dat <- matrix(abs(rnorm(25, 2, 1)), nrow = 5, ncol = 5)
class(dat) <- c("gaussian", class(dat)[length(class(dat))])
dat_org = dat

svd_res <- svd(dat)
pred_mat2 <- svd_res$u[,1:2]%*%diag(svd_res$d[1:2])%*%t(svd_res$v[,1:2])
pred_mat2

init <- initialization(dat, family = "gaussian", max_val = 100, k = 2)
pred_mat1 <- init$u_mat %*% t(init$v_mat)
pred_mat1

################

k = 2
family = "gaussian"
max_val = 100
max_iter = 10
tol = 1e-3
verbose = F

direction <- .dictate_direction(family)
if(!is.na(direction) & !is.na(max_val)){
  stopifnot((direction == ">=" & max_val > 0) | (direction == "<=" & max_val < 0))
}

# initialize
dat <- .matrix_completion(dat, k = k)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

# projected gradient descent
# pred_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
#                                         direction = direction,
#                                         max_iter = max_iter,
#                                         tol = tol)

svd_res <- svd(dat)
(svd_res$u[,1:2]%*%diag(svd_res$d[1:2])%*%t(svd_res$v[,1:2]))[1,1]


#################

n <- nrow(dat); d <- ncol(dat)
# pred_mat <- .determine_initial_matrix(dat, class(dat)[1], k = k, max_val = max_val)
# iter <- 1
# new_obj <- .evaluate_objective_mat(dat, pred_mat)
# old_obj <- Inf

#################

dat[which(dat <= tol)] <- tol/2
pred_mat <- .mean_transformation(dat, family)
direction <- .dictate_direction(family)


svd_res <- svd(dat)
(svd_res$u[,1:2]%*%diag(svd_res$d[1:2])%*%t(svd_res$v[,1:2]))[1,1]

# .project_rank_feasibility(pred_mat, k = k, direction = direction,
#                           max_val = max_val)$matrix

####################

mat <- pred_mat
stopifnot(!is.na(max_val) | !is.na(direction))
if(!is.na(max_val) & !is.na(direction)) stopifnot((direction == "<=" & max_val < 0) | (direction == ">=" & max_val > 0))

mat_org <- mat
iter <- 1

res <- .svd_projection(mat, k = k, factors = T)
res$u_mat %*% t(res$v_mat)

.svd_projection(mat, k = k, factors = F)

#######################


if(min(dim(mat)) > k+2){
  res <- tryCatch({
    # ask for more singular values than needed to ensure stability
    RSpectra::svds(mat, k = k + 2)
  }, error = function(e){
    svd(mat)
  })
} else {
  res <- svd(mat)
}

if(k == 1){
  diag_mat <- matrix(res$d[1], 1, 1)
} else {
  diag_mat <- diag(res$d[1:k])
}
res$u[,1:k,drop = F] %*% diag_mat %*% t(res$v[,1:k,drop = F])


res2 <- svd(mat)

sum(abs(res$d - res2$d[1:4]))
sum(abs(res$u - res2$u[,1:4]))
sum(abs(res$v - res2$v[,1:4]))

plot(as.numeric(res$u), as.numeric(res2$u[,1:4]), asp = T)
plot(as.numeric(res$v), as.numeric(res2$v[,1:4]), asp = T)

res3 <- res
res3$u[,2] <- -res3$u[,2]
res3$v[,2] <- -res3$v[,2]

plot(as.numeric(res3$u), as.numeric(res2$u[,1:4]), asp = T)
plot(as.numeric(res3$v), as.numeric(res2$v[,1:4]), asp = T)

# res$u[,1:k,drop = F] %*% diag_mat %*% t(res$v[,1:k,drop = F])
res2$u[,1:k,drop = F] %*% diag(res2$d[1:k]) %*% t(res2$v[,1:k,drop = F])



