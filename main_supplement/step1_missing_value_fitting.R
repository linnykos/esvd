set.seed(10)
load(paste0("../results/baron_step0_screening", suffix, ".RData"))

neg_binom_vec <- c(50, 100, 500, 1000, 5000, 1e4, 1e6)
curved_gaussian_vec <- c(0.1, 0.5, 1, 1.5, 2, 4, 100)

max_iter <- 50
max_val <- 3000
paramMat <- as.matrix(expand.grid(1:4, 3:5, 50, 3000, 10))
colnames(paramMat) <- c("fitting_distr", "k")

######################################

fitting_func <- function(dat_impute, fitting_distr, k, missing_idx){
  dat_NA <- dat_impute
  dat_NA[missing_idx] <- NA
  set.seed(10)

  if(fitting_distr %in% c(1,2)){
    fit_list <- vector("list", 1)

    set.seed(10)
    init <- eSVD::initialization(dat_NA, family = family, k = k, max_val = max_val)
    fit_list[[1]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                        family = family,
                                        max_iter = max_iter, max_val = max_val,
                                        return_path = F, cores = ncores,
                                        verbose = F)
  } else if(fit == 3) {
    fit_list <- lapply(neg_binom_vec, function(scalar){
      set.seed(10)
      init <- eSVD::initialization(dat_NA, family = "neg_binom", k = k, max_val = max_val,
                                   scalar = scalar)
      fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                     family = "neg_binom", scalar = scalar,
                                     max_iter = max_iter, max_val = max_val,
                                     return_path = F, cores = ncores,
                                     verbose = F)

      fit
    })
  } else {
    fit_list <- lapply(curved_gaussian_vec, function(scalar){
      set.seed(10)
      init <- eSVD::initialization(dat_NA, family = "curved_gaussian", k = k, max_val = max_val,
                                   scalar = scalar)
      fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                     family = "curved_gaussian", scalar = scalar,
                                     max_iter = max_iter, max_val = max_val,
                                     return_path = F, cores = ncores,
                                     verbose = F)

      fit
    })
  }
}

res_missing_list <- vector("list", length(preprocessing_list))
for(i in 1:length(preprocessing_list)){
  print(paste0("On dataset ", i))
  set.seed(10)
  missing_idx <- eSVD::construct_missing_values(n = nrow(preprocessing_list[[i]]$dat_impute),
                                                p = ncol(preprocessing_list[[i]]$dat_impute),
                                                num_val = 2)

  tmp_res <- lapply(1:nrow(paramMat), function(j){
    print(paste0("On row ", j))
    fitting_func(preprocessing_list[[i]]$dat_impute, fitting_distr = paramMat[j, "fitting_distr"],
                 k = paramMat[j, "k"], missing_idx = missing_idx)
  })

  res_missing_list[[i]] <- list(fit_list = tmp_res, missing_idx = missing_idx)
}

rm(list = c("tmp_res", "missing_idx"))
save.image(paste0("../results/baron_step1_missing_value_fitting", suffix, ".RData"))

