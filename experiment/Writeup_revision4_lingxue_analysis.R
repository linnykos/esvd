rm(list=ls())

labels_file_vec <- c("/raid6/Kevin/singlecell_data/Cleaned/Baron_human1/Baron_human1_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_human2/Baron_human2_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_human3/Baron_human3_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_human4/Baron_human4_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse1/Baron_mouse1_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse2/Baron_mouse2_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Darmanis/Darmanis_label.csv")
dat_file_vec <- c("/raid6/Kevin/singlecell_data/Cleaned/Baron_human1/Baron_human1_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_human2/Baron_human2_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_human3/Baron_human3_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_human4/Baron_human4_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse1/Baron_mouse1_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse2/Baron_mouse2_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Darmanis/Darmanis_expr.csv")

ncores <- 20
doMC::registerDoMC(cores = ncores)

preprocessing_list <- vector("list", length(dat_file_vec))
for(i in 1:length(labels_file_vec)){
  labels <- read.csv(labels_file_vec[i])
  dat <- read.csv(dat_file_vec[i])
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  label_vec <- labels$cluster

  ##########

  tab <- table(label_vec)
  idx <- which(tab/sum(tab) <= 0.05)

  dat <- dat[which(!label_vec %in% names(tab[idx])),]
  label_vec <- label_vec[which(!label_vec %in% names(tab[idx]))]
  label_vec <- as.factor(as.character(label_vec))

  # next remove genes
  tmp <- apply(dat, 2, function(x){sum(x != 0)})
  dat <- dat[,which(tmp >= floor(nrow(dat)/100))]
  dat <- as.matrix(dat)

  # try a series of SPCAs
  k <- 5
  lvls <- 20
  v_seq <- exp(seq(log(1), log(sqrt(ncol(dat))/2), length.out = lvls))
  res_list <- vector("list", lvls)

  spca_func <- function(i){
    res <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
    print(paste0("Finished SPC for level ", i))
    res
  }

  res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))

  spca_mat <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
  idx <- min(which(spca_mat[,2] == ncol(dat)))
  target_var <- spca_mat[idx,3]
  idx <- min(intersect(which(spca_mat[,2] >= 15), which(spca_mat[,3] >= 0.8*target_var)))
  spca_idx <- sort(unlist(apply(res_list[[idx]]$v, 2, function(x){which(x != 0)})))

  res_descend <- descend::runDescend(t(dat), n.cores = ncores)
  res_hvg <- descend::findHVG(res_descend, threshold = 50)
  descend_idx <- which(colnames(dat) %in% res_hvg$HVG.genes)

  gene_idx <- sort(unique(c(spca_idx, descend_idx)))
  dat_impute <- dat[,gene_idx]

  reweight_factor <- rowSums(dat_impute)
  dat_impute <- t(sapply(1:nrow(dat_impute), function(x){dat_impute[x,]/reweight_factor[x]}))
  dat_impute <- dat_impute * 1000/max(dat_impute)

  preprocessing_list[[i]] <- list(dat_impute = dat_impute, label_vec = label_vec,
                                  gene_idx = gene_idx, spca_mat = spca_mat)

  save.image("../results/lingxue_analysis.RData")
}

rm(list=c("dat"))
save.image("../results/lingxue_analysis.RData")

fit_all_list <- vector("list", length(preprocessing_list))
## now do all the fits
for(num in 1:length(preprocessing_list)){
  dat_impute <- preprocessing_list[[num]]$dat_impute

  fit_list <- vector("list", 12)
  names(fit_list) <- c("gaussian_missing", "gaussian", "poisson_missing", "poisson",
                       "exponential_missing", "exponential", "neg_binom_missing",
                       "neg_binom", "curved_gaussian_missing", "curved_gaussian",
                       "neg_bin_param", "curved_gaussian_param")

  n <- nrow(dat_impute); d <- ncol(dat_impute)
  set.seed(10)
  missing_idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
    cbind(x, sample(1:d, 2))
  }))), do.call(rbind, (lapply(1:d, function(x){
    cbind(sample(1:n, 2), x)
  }))))

  dat_NA <- dat_impute
  for(i in 1:nrow(missing_idx)){
    dat_NA[missing_idx[i,1], missing_idx[i,2]] <- NA
  }
  missing_idx <- which(is.na(dat_NA))
  missing_val <- dat_impute[missing_idx]

  # parameters to fit
  max_iter <- 50
  max_val <- 3000
  k <- 5
  fitting_iter <- 5

  # gaussian fit
  init <- eSVD::initialization(dat_NA, family = "gaussian", k = k, max_val = max_val)
  fit_list[[1]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                  family = "gaussian",
                                                  max_iter = max_iter, max_val = max_val,
                                                  return_path = F, cores = ncores,
                                                  verbose = F)
  save.image("../results/lingxue_analysis.RData")

  init <- eSVD::initialization(dat_impute, family = "gaussian", k = k, max_val = max_val)
  fit_list[[2]] <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                          family = "gaussian",
                                          max_iter = max_iter, max_val = max_val,
                                          return_path = F, cores = ncores,
                                          verbose = F)
  save.image("../results/lingxue_analysis.RData")

  # poisson fit
  init <- eSVD::initialization(dat_NA, family = "poisson", k = k, max_val = max_val)
  fit_list[[3]]  <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = "poisson",
                                                 max_iter = max_iter, max_val = max_val,
                                                 return_path = F, cores = ncores,
                                                 verbose = F)
  save.image("../results/lingxue_analysis.RData")

  init <- eSVD::initialization(dat_impute, family = "poisson", k = k, max_val = max_val)
  fit_list[[4]] <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                         family = "poisson",
                                         max_iter = max_iter, max_val = max_val,
                                         return_path = F, cores = ncores,
                                         verbose = F)
  save.image("../results/lingxue_analysis.RData")

  # exponential fit
  init <- eSVD::initialization(dat_NA, family = "exponential", k = k, max_val = max_val)
  fit_list[[5]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                     family = "exponential",
                                                     max_iter = max_iter, max_val = max_val,
                                                     return_path = F, cores = ncores,
                                                     verbose = F)
  save.image("../results/lingxue_analysis.RData")

  init <- eSVD::initialization(dat_impute, family = "exponential", k = k, max_val = max_val)
  fit_list[[6]] <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                             family = "exponential",
                                             max_iter = max_iter, max_val = max_val,
                                             return_path = F, cores = ncores,
                                             verbose = F)
  save.image("../results/lingxue_analysis.RData")

  # negative binomial fit
  neg_bin_param <- eSVD::tuning(dat_impute, fit_list[[4]]$u_mat, fit_list[[4]]$v_mat, family_to = "neg_binom",
                                family_from = "poisson")
  for(i in 1:fitting_iter){
    init <- eSVD::initialization(dat_impute, family = "neg_binom", k = k, max_val = max_val,
                                 size = neg_bin_param)
    fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "neg_binom", size = neg_bin_param,
                                   max_iter = max_iter, max_val = max_val,
                                   return_path = F, cores = ncores,
                                   verbose = F)

    neg_bin_param <- eSVD::tuning(dat_impute, fit$u_mat, fit$v_mat, family_to = "neg_binom",
                                  family_from = "neg_binom")
  }
  save.image("../results/lingxue_analysis.RData")

  init <- eSVD::initialization(dat_NA, family = "neg_binom", k = k, max_val = max_val,
                               size = neg_bin_param)
  fit_list[[7]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = "neg_binom", size = neg_bin_param,
                                                   max_iter = max_iter, max_val = max_val,
                                                   return_path = F, cores = ncores,
                                                   verbose = F)
  save.image("../results/lingxue_analysis.RData")

  init <- eSVD::initialization(dat_impute, family = "neg_binom", k = k, max_val = max_val,
                               size = neg_bin_param)
  fit_list[[8]] <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                           family = "neg_binom", size = neg_bin_param,
                                           max_iter = max_iter, max_val = max_val,
                                           return_path = F, cores = ncores,
                                           verbose = F)
  save.image("../results/lingxue_analysis.RData")

  # curved gaussian fit
  curved_gaussian_param <- eSVD::tuning(dat_impute, fit_list[[6]]$u_mat, fit_list[[6]]$v_mat, family_to = "curved_gaussian",
                                        family_from = "exponential")
  for(i in 1:fitting_iter){
    init <- eSVD::initialization(dat_impute, family = "curved_gaussian", k = k, max_val = max_val,
                                 scalar = curved_gaussian_param)
    fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "curved_gaussian", scalar = curved_gaussian_param,
                                   max_iter = max_iter, max_val = max_val,
                                   return_path = F, cores = ncores,
                                   verbose = F)


    curved_gaussian_param <- eSVD::tuning(dat_impute, fit$u_mat, fit$v_mat, family_to = "curved_gaussian",
                                          family_from = "curved_gaussian")
  }
  save.image("../results/lingxue_analysis.RData")

  init <- eSVD::initialization(dat_NA, family = "curved_gaussian", k = k, max_val = max_val,
                               scalar = curved_gaussian_param)
  fit_list[[9]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                         family = "curved_gaussian", scalar = curved_gaussian_param,
                                                         max_iter = max_iter, max_val = max_val,
                                                         return_path = F, cores = ncores,
                                                         verbose = F)
  save.image("../results/lingxue_analysis.RData")

  init <- eSVD::initialization(dat_impute, family = "curved_gaussian", k = k, max_val = max_val,
                               scalar = curved_gaussian_param)
  fit_list[[10]] <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = "curved_gaussian", scalar = curved_gaussian_param,
                                                 max_iter = max_iter, max_val = max_val,
                                                 return_path = F, cores = ncores,
                                                 verbose = F)

  fit_list[[11]] <- neg_bin_param
  fit_list[[12]] <-curved_gaussian_param

  fit_all_list[[num]] <- fit_list

  save.image("../results/lingxue_analysis.RData")
}

save.image("../results/lingxue_analysis.RData")
