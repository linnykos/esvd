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

# preprocessing_list <- vector("list", length(dat_file_vec))
# for(i in 1:length(labels_file_vec)){
#   labels <- read.csv(labels_file_vec[i])
#   dat <- read.csv(dat_file_vec[i])
#   dat <- dat[,-1]
#   dat <- as.matrix(dat)
#   label_vec <- labels$cluster
#
#   ##########
#
#   tab <- table(label_vec)
#   idx <- which(tab/sum(tab) <= 0.05)
#
#   dat <- dat[which(!label_vec %in% names(tab[idx])),]
#   label_vec <- label_vec[which(!label_vec %in% names(tab[idx]))]
#   label_vec <- as.factor(as.character(label_vec))
#
#   # next remove genes
#   tmp <- apply(dat, 2, function(x){sum(x != 0)})
#   dat <- dat[,which(tmp >= floor(nrow(dat)/100))]
#   dat <- as.matrix(dat)
#
#   # try a series of SPCAs
#   k <- 5
#   lvls <- 20
#   v_seq <- exp(seq(log(1), log(sqrt(ncol(dat))/2), length.out = lvls))
#   res_list <- vector("list", lvls)
#
#   spca_func <- function(i){
#     res <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
#     print(paste0("Finished SPC for level ", i))
#     res
#   }
#
#   res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))
#
#   spca_mat <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
#   idx <- min(which(spca_mat[,2] == ncol(dat)))
#   target_var <- spca_mat[idx,3]
#   idx <- min(intersect(which(spca_mat[,2] >= 15), which(spca_mat[,3] >= 0.8*target_var)))
#   spca_idx <- sort(unlist(apply(res_list[[idx]]$v, 2, function(x){which(x != 0)})))
#
#   res_descend <- descend::runDescend(t(dat), n.cores = ncores)
#   res_hvg <- descend::findHVG(res_descend, threshold = 50)
#   descend_idx <- which(colnames(dat) %in% res_hvg$HVG.genes)
#
#   gene_idx <- sort(unique(c(spca_idx, descend_idx)))
#   dat_impute <- dat[,gene_idx]
#
#   reweight_factor <- rowSums(dat_impute)
#   dat_impute <- t(sapply(1:nrow(dat_impute), function(x){dat_impute[x,]/reweight_factor[x]}))
#   dat_impute <- dat_impute * 1000/max(dat_impute)
#
#   preprocessing_list[[i]] <- list(dat_impute = dat_impute, label_vec = label_vec,
#                                   gene_idx = gene_idx, spca_mat = spca_mat)
#
#   save.image("../results/lingxue_data_preprocessed.RData")
# }
#
# rm(list=c("dat"))
# save.image("../results/lingxue_data_preprocessed.RData")
load("../results/lingxue_data_preprocessed.RData")

sapply(preprocessing_list, function(x){dim(x$dat_impute)})

######################################

fitting_func <- function(dat_impute, k, missing_idx){
  max_iter <- 50
  max_val <- 3000
  iter_max <- 10

  dat_NA <- dat_impute
  dat_NA[missing_idx] <- NA

  fit_list <- vector("list", 8)
  names(fit_list) <- c("gaussian_missing", "poisson_missing", "neg_binom_missing",
                       "curved_gaussian_missing", "neg_bin_param", "curved_gaussian_param",
                       "missing_idx", "missing_val")

  # gaussian fit
  print("Gaussian")
  init <- eSVD::initialization(dat_NA, family = "gaussian", k = k, max_val = max_val)
  fit_list[[1]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                           family = "gaussian",
                                           max_iter = max_iter, max_val = max_val,
                                           return_path = F, cores = ncores,
                                           verbose = F)
  save.image("../results/lingxue_analysis_tmp.RData")

  #################

  # poisson fit
  print("Poisson")
  init <- eSVD::initialization(dat_NA, family = "poisson", k = k, max_val = max_val)
  fit_list[[2]]  <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                            family = "poisson",
                                            max_iter = max_iter, max_val = max_val,
                                            return_path = F, cores = ncores,
                                            verbose = F)
  save.image("../results/lingxue_analysis_tmp.RData")

  #################

  # negative binomial fit
  print("Neg binomial")
  neg_bin_vec <- eSVD::tuning_scalar(dat_impute, family = "neg_binom",
                                     max_iter = max_iter, max_val = max_val, k = k,
                                     return_path = F, cores = ncores, iter_max = iter_max,
                                     search_min = 1, search_max = 4000)
  neg_bin_param <- neg_bin_vec[length(neg_bin_vec)]
  save.image("../results/lingxue_analysis_tmp.RData")

  init <- eSVD::initialization(dat_NA, family = "neg_binom", k = k, max_val = max_val,
                               scalar = neg_bin_param)
  fit_list[[3]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                           family = "neg_binom", scalar = neg_bin_param,
                                           max_iter = max_iter, max_val = max_val,
                                           return_path = F, cores = ncores,
                                           verbose = F)
  save.image("../results/lingxue_analysis_tmp.RData")

  #################

  # curved gaussian fit
  print("Curved gaussian")
  curved_gaussian_vec <- eSVD::tuning_scalar(dat_impute, family = "curved_gaussian",
                                             max_iter = max_iter, max_val = max_val, k = k,
                                             return_path = F, cores = ncores, iter_max = iter_max,
                                             search_min = 0.1, search_max = 10)
  curved_gaussian_param <- curved_gaussian_vec[length(curved_gaussian_vec)]
  save.image("../results/lingxue_analysis_tmp.RData")

  init <- eSVD::initialization(dat_NA, family = "curved_gaussian", k = k, max_val = max_val,
                               scalar = curved_gaussian_param)
  fit_list[[4]] <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                           family = "curved_gaussian", scalar = curved_gaussian_param,
                                           max_iter = max_iter, max_val = max_val,
                                           return_path = F, cores = ncores,
                                           verbose = F)
  save.image("../results/lingxue_analysis_tmp.RData")

  fit_list[[5]] <- neg_bin_vec
  fit_list[[6]] <- curved_gaussian_vec
  fit_list[[7]] <- missing_idx
  fit_list[[8]] <- dat_impute[missing_idx]

  fit_list
}

######

fit_all_list <- vector("list", length(preprocessing_list))

## now do all the fits
for(num in 1:length(preprocessing_list)){
  print(paste0("Working on dataset ", num))
  dat_impute <- preprocessing_list[[num]]$dat_impute

  set.seed(10)
  missing_idx <- eSVD::construct_missing_values(n = nrow(dat_impute), p = ncol(dat_impute), num_val = 2)

  # parameters to fit
  k_vec <- 3:5

  fit_list <- vector("list", length(k_vec))
  for(i in 1:3){
    print(paste0("Working on dimension ", k_vec[i]))
    fit_list[[i]] <-  fitting_func(dat_impute, k_vec[i], missing_idx)
    save.image("../results/lingxue_analysis.RData")
  }

  fit_all_list[[num]] <- fit_list
  save.image("../results/lingxue_analysis.RData")
}

save.image("../results/lingxue_analysis.RData")
