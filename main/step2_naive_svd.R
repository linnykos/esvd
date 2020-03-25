set.seed(10)
load(paste0("../results/step1_imputing", suffix, ".RData"))

max_val <- 5000
n <- nrow(dat_impute); d <- ncol(dat_impute)
cv_trials <- 3
svd_missing_list <- vector("list", length(cv_trials))
k <- 5

# generate list of missing indices
missing_idx_list <- lapply(1:cv_trials, function(j){
  set.seed(10*j)
  eSVD::construct_missing_values(n = n, p = d, num_val = 4)
})

# fit many softImputes
starting_lambda <- min(sapply(1:cv_trials, function(i){
  dat_NA <- dat_impute
  dat_NA[missing_idx_list[[i]]] <- NA

  softImpute::lambda0(dat_NA)
}))

paramMat_svd <- as.matrix(expand.grid(seq(1, starting_lambda, length.out = 50), c(3,4,5)))
colnames(paramMat_svd) <- c("lambda", "k")

svd_angle_vec <- sapply(1:nrow(paramMat_svd), function(i){
  tmp_mat <- do.call(rbind, lapply(1:cv_trials, function(i){
    dat_NA <- dat_impute
    dat_NA[missing_idx_list[[i]]] <- NA

    softImpute_embedding <- softImpute::softImpute(dat_NA, rank.max = paramMat[i, "k"], lambda = paramMat[i, "lambda"])
    softImpute_pred <- softImpute_embedding$u %*% diag(softImpute_embedding$d) %*% t(softImpute_embedding$v)
    softImpute_pred <- softImpute_pred[missing_idx_list[[i]]]
    obs_val <- dat_impute[missing_idx_list[[i]]]

    cbind(obs_val, softImpute_pred)
  }))

  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))
  ang * 180/pi
})

# select the best model
idx <- which.min(abs(svd_angle_vec - 45))
svd_missing <- lapply(1:cv_trials , function(i){
  softImpute::softImpute(dat_NA, rank.max = paramMat[idx, "k"], lambda = paramMat[idx, "lambda"])
})

# fit the actual embedding as well
tmp <- svd(dat_impute)
svd_embedding <- tmp$u[,1:paramMat[idx, "k"]] %*% diag(sqrt(tmp$d[1:paramMat[idx, "k"]]))

rm(list = c("init", "fit", "dat_impute_NA", "j", "tmp", "starting_lambda", "idx"))
print(paste0(Sys.time(), ": Finished naive SVD"))
save.image(paste0("../results/step2_naive_svd", suffix, ".RData"))
