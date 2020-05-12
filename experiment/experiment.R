rm(list=ls())
load("../results/tmp2.RData")

i <- 3
ii <- 6
dat_impute = preprocessing_list[[i]]$dat_impute
vec = paramMat_esvd[ii,]
missing_idx_list = missing_idx_list_list[[i]]

fitting_distr <- c("gaussian", "neg_binom", "curved_gaussian")[vec["fitting_distr"]]

dat_org <- dat_impute*1000/max(dat_impute)
# set missing value
dat_NA <- dat_org

for(jj in 1:3){
  dat_NA[missing_idx_list[[jj]]] <- NA
}

set.seed(10)
init <- eSVD::initialization(dat_NA, family = fitting_distr, k = vec["k"], max_val = vec["max_val"],
                             scalar = vec["scalar"])
# zz <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
#                                                 family = fitting_distr, max_iter = vec["max_iter"],
#                                                 scalar = vec["scalar"],
#                                                 max_val = vec["max_val"], return_path = F, cores = ncores, verbose = T)

######################################
set.seed(10)
dat = dat_NA
k = vec["k"]
max_val = vec["max_val"]
scalar = vec["scalar"]
family = fitting_distr
max_val = NA
max_iter = 10
tol = 1e-3
verbose = F

direction <- eSVD:::.dictate_direction(family)
if(!is.na(max_val)){
  if(!is.na(direction) && direction == "<=") max_val <- -max_val
}

# initialize
dat <- eSVD:::.matrix_completion(dat, k = k)
dat[1:5,1:5]
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

# projected gradient descent
nat_mat <- eSVD:::.projected_gradient_descent(dat, k = k, max_val = max_val,
                                       direction = direction,
                                       max_iter = max_iter,
                                       tol = tol, scalar = scalar)
nat_mat[1:5,1:5]

res <- eSVD:::.svd_projection(nat_mat, k = k, factors = T)
u_mat <- res$u_mat; v_mat <- res$v_mat
head(u_mat)

if(!is.na(direction)){
  if(direction == "<=") {
    stopifnot(all(nat_mat[which(!is.na(dat))] < 0))
  } else {
    stopifnot(all(nat_mat[which(!is.na(dat))] > 0))
  }
}

tmp <- eSVD:::.fix_rank_defficiency_initialization(u_mat, v_mat, direction)
u_mat <- tmp$u_mat; v_mat <- tmp$v_mat
head(u_mat)

# tmp <- eSVD:::.reparameterize(u_mat, v_mat)
# u_mat <- tmp$u_mat; v_mat <- tmp$v_mat
# head(u_mat)

#########################
check = F
tol = 1e-6
cov_x <- t(u_mat) %*% u_mat
cov_y <- t(v_mat) %*% v_mat

stopifnot(all(dim(cov_x) == dim(cov_y)), nrow(cov_x) == ncol(cov_x))
if(nrow(cov_x) == 1){
  return(matrix((as.numeric(cov_y)/as.numeric(cov_x))^(1/4), 1, 1))
}

eigen_x <- eigen(cov_x)
eigen_y <- eigen(cov_y)

Vx <- eigen_x$vectors
Vy <- eigen_y$vectors

if(any(eigen_x$values <= tol) | any(eigen_y$values <= tol)) warning("Detecting rank defficiency in reparameterization step")

Dx <- diag(eigen_x$values)
Dy <- diag(eigen_y$values)

tmp <- sqrt(Dy) %*% t(Vy) %*% Vx %*% sqrt(Dx)
svd_tmp <- svd(tmp)
R <- svd_tmp$u %*% t(svd_tmp$v)

Dx_inv <- Dx; diag(Dx_inv) <- 1/diag(Dx)
sym_prod <- Vx %*% sqrt(Dx_inv) %*% t(R) %*% sqrt(Dy) %*% t(Vy)
sym_prod[which(abs(sym_prod) <= tol)] <- 0

eigen_sym <- eigen(sym_prod)
eigen_sym
