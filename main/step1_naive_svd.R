set.seed(10)
load(paste0("../results/step0_screening", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

max_val <- 5000
n <- nrow(dat); p <- ncol(dat)
cv_trials <- 3

# generate list of missing indices
missing_idx_list <- lapply(1:cv_trials, function(j){
  set.seed(10*j)
  eSVD::construct_missing_values(n = n, p = p, num_val = 4)
})

training_idx_list <- lapply(1:cv_trials, function(j){
  c(1:prod(dim(dat)))[-missing_idx_list[[j]]]
})

log_dat <- log2(dat+1)

# fit many softImputes
starting_lambda <- min(sapply(1:cv_trials, function(i){
  log_dat_NA <- log_dat
  log_dat_NA[missing_idx_list[[i]]] <- NA

  softImpute::lambda0(log_dat_NA)
}))

paramMat_svd <- as.matrix(expand.grid(seq(1, starting_lambda, length.out = 50), c(5,10,20,50)))
colnames(paramMat_svd) <- c("lambda", "k")

svd_angle_res <- sapply(1:nrow(paramMat_svd), function(j){
  tmp_res <- sapply(1:cv_trials, function(i){
    log_dat_NA <- log_dat
    log_dat_NA[missing_idx_list[[i]]] <- NA

    set.seed(10)
    softImpute_embedding <- softImpute::softImpute(log_dat_NA, rank.max = paramMat_svd[j, "k"], lambda = paramMat_svd[j, "lambda"])
    softImpute_pred <- softImpute_embedding$u %*% diag(softImpute_embedding$d) %*% t(softImpute_embedding$v)

    tmp_mat <- cbind(as.numeric(log_dat), as.numeric(softImpute_pred))

    training_val <- eSVD::compute_principal_angle(tmp_mat[training_idx_list[[i]],])
    testing_val <- eSVD::compute_principal_angle(tmp_mat[missing_idx_list[[i]],])

    c(training = training_val, testing = testing_val)
  })

  rowMeans(tmp_res)
})
svd_angle_res <- cbind(t(svd_angle_res), paramMat_svd)

# select the best model
idx <- which.min(abs(svd_angle_res[,2] - 45))
svd_missing_list <- lapply(1:cv_trials , function(i){
  log_dat_NA <- log_dat
  log_dat_NA[missing_idx_list[[i]]] <- NA

  set.seed(10)
  softImpute::softImpute(log_dat_NA, rank.max = paramMat_svd[idx, "k"], lambda = paramMat_svd[idx, "lambda"])
})


# fit the actual embedding as well
tmp <- svd(log_dat)
svd_scree <- tmp$d
svd_embedding <- (n/p)^(1/4)*tmp$u[,1:paramMat_svd[idx, "k"]] %*% diag(sqrt(tmp$d[1:paramMat_svd[idx, "k"]]))

rm(list = c("init", "fit", "dat_NA", "j", "tmp", "starting_lambda", "idx",
            "log_dat", "log_dat_NA"))
print(paste0(Sys.time(), ": Finished naive SVD"))
source_code_info <- c(source_code_info, readLines("../main/step1_naive_svd.R"))
save.image(paste0("../results/step1_naive_svd", suffix, ".RData"))
print(warnings())
