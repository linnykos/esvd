rm(list=ls())
zz <- .mean_transformation(dat_impute, family = "curved_gaussian", scalar = scalar_val)
svd_res <- svd(zz)
plot(svd_res$d)

# zz <- eigen(cor(dat_impute))
# plot(zz$values)

a <- .evaluate_objective_mat.curved_gaussian(dat_impute, res_our$u_mat %*% t(res_our$v_mat), scalar = scalar_val)
b <- .evaluate_objective_mat.curved_gaussian(dat_impute, zz, scalar = scalar_val)
tmp_vec <- rowMeans(zz)
tmp_mat <- tmp_vec %*% t(rep(1, ncol(dat_impute)))
c <- .evaluate_objective_mat.curved_gaussian(dat_impute, tmp_mat, scalar = scalar_val)

(a-c)/(b-c)
