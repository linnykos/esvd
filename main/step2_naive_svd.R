set.seed(10)
load(paste0("../results/step1_imputing", suffix, ".RData"))

max_val <- 5000
n <- nrow(dat_impute); d <- ncol(dat_impute)
cv_trials <- 5
naive_res_list <- vector("list", length(cv_trials))
k <- 5

# generate list of missing indices
missing_idx_list <- lapply(1:cv_trials, function(j){
  set.seed(10*j)
  eSVD::construct_missing_values(n = n, p = d, num_val = 4)
})

for(j in 1:cv_trials){
  print(paste0("On trial ", j))

  # set missing values
  dat_impute_NA <- dat_impute
  dat_impute_NA[missing_idx_list[[j]]] <- NA

  # fit
  set.seed(10)
  init <- eSVD::initialization(dat_impute_NA, family = "gaussian",
                               k = k, max_val = max_val)
  fit <- eSVD::fit_factorization(dat_impute_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                family = "gaussian",
                                                max_iter = 50, max_val = max_val,
                                                return_path = F, cores = ncores,
                                                verbose = T)

  pred_mat <- fit$u_mat %*% t(fit$v_mat)
  pred_val <- pred_mat[missing_idx_list[[j]]]
  observed_val <- dat_impute[missing_idx_list[[j]]]

  naive_res_list[[j]] <- cbind(pred_val, observed_val)

  save.image(paste0("../results/step2_naive_svd", suffix, "_tmp.RData"))
}

tmp <- svd(dat_impute)
res_svd <- tmp$u[,1:k] %*% diag(sqrt(tmp$d[1:k]))

rm(list = c("init", "fit", "dat_impute_NA", "j", "tmp"))
print(paste0(Sys.time(), ": Finished naive SVD"))
save.image(paste0("../results/step2_naive_svd", suffix, ".RData"))
