var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step7_figures", suffix, ".RData"))

nat_mat_list <- lapply(1:length(svd_missing_list), function(i){
  svd_missing_list[[i]]$u %*% diag(svd_missing_list[[i]]$d) %*% t(svd_missing_list[[i]]$v)
})

tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
  cbind(log2(dat_impute/rescaling_factor + 1)[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]])
}))
sd_val <- sd(tmp_mat[,1] - tmp_mat[,2])

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

png(filename = paste0("../../esvd_results/figure/main/svd_training_testing.png"),
    height = 1500, width = 2500, res = 300,
    units = "px")
par(mfrow = c(1,2))
eSVD::plot_prediction_against_observed(log2(dat_impute/rescaling_factor+1), nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "gaussian", scalar = sd_val,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 transparency = 0.05, xlim = c(-2,15), ylim = c(-2,15),
                                 max_points = 1e6, cex.lab = 1.25)

eSVD::plot_prediction_against_observed(log2(dat_impute/rescaling_factor+1), nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "gaussian", scalar = sd_val,  transparency = 0.05,
                                 xlim = c(-2,15), ylim = c(-2,15),
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Testing set)",
                                 cex.lab = 1.25)
graphics.off()

############################

nat_mat_list <- lapply(1:length(esvd_missing_list[[esvd_angle_res$idx]]), function(i){
  esvd_missing_list[[esvd_angle_res$idx]][[i]]$u_mat %*% t(esvd_missing_list[[esvd_angle_res$idx]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

png(filename = paste0("../../esvd_results/figure/main/esvd_training_testing.png"),
    height = 1500, width = 2500, res = 300,
    units = "px")
par(mfrow = c(1,2))
eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 transparency = 0.2, cex.lab = 1.25,
                                 family = "curved_gaussian",
                                 scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 transparency = 0.2, cex.lab = 1.25,
                                 family = "curved_gaussian",
                                 scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")
graphics.off()

png(filename = paste0("../../esvd_results/figure/main/esvd_training_testing_zoomedin.png"),
    height = 1500, width = 2500, res = 300,
    units = "px")
par(mfrow = c(1,2))
eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                       missing_idx_list = training_idx_list,
                                       family = "curved_gaussian", xlim = c(0,50), ylim = c(0,50),
                                       transparency = 0.05,
                                       scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                       main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                       max_points = 1e6)


eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                       missing_idx_list = missing_idx_list,
                                       family = "curved_gaussian", xlim = c(0,50), ylim = c(0,50),
                                       transparency = 0.05,
                                       scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                       main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")
graphics.off()

