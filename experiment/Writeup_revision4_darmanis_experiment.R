rm(list=ls())
labels_file <- "/raid6/Kevin/singlecell_data/Cleaned/Darmanis/Darmanis_label.csv"
dat_file <- "/raid6/Kevin/singlecell_data/Cleaned/Darmanis/Darmanis_expr.csv"

ncores <- 10
doMC::registerDoMC(cores = ncores)

########

labels <- read.csv(labels_file)
dat <- read.csv(dat_file)
dat <- dat[,-1]
dat <- as.matrix(dat)
label_vec <- labels$cluster

##########

tab <- table(label_vec)
idx <- which(tab/sum(tab) <= 0.05)

dat <- dat[which(!label_vec %in% names(tab[idx])),]
label_vec <- label_vec[which(!label_vec %in% names(tab[idx]))]

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

print("Starting SPCA")
res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))

tmp_spca_mat <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
idx_spca <- which.min(abs(tmp_spca_mat[,3] - 0.9))
if(length(idx_spca) == 0) idx_spca <- nrow(tmp_spca_mat)

print("Starting DESCEND")
res_descend <- descend::runDescend(t(dat), n.cores = ncores)
res_hvg <- descend::findHVG(res_descend, threshold = 50)

idx1 <- sort(unlist(apply(res_list[[idx_spca]]$v, 2, function(x){which(x != 0)})))
idx2 <- which(colnames(dat) %in% res_hvg$HVG.genes)
idx <- sort(unique(c(idx1, idx2)))

dat_impute <- dat[,idx]

save.image("../results/darmanis.RData")

# reweight rows
reweight_factor <- rowSums(dat_impute)
dat_impute <- t(sapply(1:nrow(dat_impute), function(x){dat_impute[x,]/reweight_factor[x]}))
dat_impute <- dat_impute * 1000/max(dat_impute)

######################

# make missing values
n <- nrow(dat_impute); d <- ncol(dat_impute)
missing_idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
  cbind(x, sample(1:d, 2))
}))), do.call(rbind, (lapply(1:d, function(x){
  cbind(sample(1:n, 2), d)
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
gaussian_fit_missing <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "gaussian",
                               max_iter = max_iter, max_val = max_val,
                               return_path = F, cores = ncores,
                               verbose = F)
save.image("../results/darmanis.RData")

init <- eSVD::initialization(dat_impute, family = "gaussian", k = k, max_val = max_val)
gaussian_fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                family = "gaussian",
                                                max_iter = max_iter, max_val = max_val,
                                                return_path = F, cores = ncores,
                                                verbose = F)
save.image("../results/darmanis.RData")

# poisson fit
init <- eSVD::initialization(dat_NA, family = "poisson", k = k, max_val = max_val)
poisson_fit_missing <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "poisson",
                               max_iter = max_iter, max_val = max_val,
                               return_path = F, cores = ncores,
                               verbose = F)
save.image("../results/darmanis.RData")

init <- eSVD::initialization(dat_impute, family = "poisson", k = k, max_val = max_val)
poisson_fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                               family = "poisson",
                                               max_iter = max_iter, max_val = max_val,
                                               return_path = F, cores = ncores,
                                               verbose = F)
save.image("../results/darmanis.RData")

# exponential fit
init <- eSVD::initialization(dat_impute, family = "exponential", k = k, max_val = max_val)
exponential_fit_missing <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "exponential",
                               max_iter = max_iter, max_val = max_val,
                               return_path = F, cores = ncores,
                               verbose = F)
save.image("../results/darmanis.RData")

init <- eSVD::initialization(dat_impute, family = "exponential", k = k, max_val = max_val)
exponential_fit <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                                                   family = "exponential",
                                                   max_iter = max_iter, max_val = max_val,
                                                   return_path = F, cores = ncores,
                                                   verbose = F)
save.image("../results/darmanis.RData")

# negative binomial fit
neg_bin_param <- eSVD::tuning(dat_impute, exponential_fit$u_mat, exponential_fit$v_mat, family_to = "neg_binom",
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
save.image("../results/darmanis.RData")

init <- eSVD::initialization(dat_NA, family = "neg_binom", k = k, max_val = max_val,
                             size = neg_bin_param)
neg_binom_fit_missing <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "neg_binom", size = neg_bin_param,
                               max_iter = max_iter, max_val = max_val,
                               return_path = F, cores = ncores,
                               verbose = F)
save.image("../results/darmanis.RData")

init <- eSVD::initialization(dat_impute, family = "neg_binom", k = k, max_val = max_val,
                             size = neg_bin_param)
neg_binom_fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                 family = "neg_binom", size = neg_bin_param,
                                                 max_iter = max_iter, max_val = max_val,
                                                 return_path = F, cores = ncores,
                                                 verbose = F)
save.image("../results/darmanis.RData")

# curved gaussian fit
curved_gaussian_param <- eSVD::tuning(dat_impute, fit$u_mat, fit$v_mat, family_to = "curved_gaussian",
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
save.image("../results/darmanis.RData")

init <- eSVD::initialization(dat_NA, family = "curved_gaussian", k = k, max_val = max_val,
                             scalar = curved_gaussian_param)
curved_gaussian_fit_missing <- eSVD::fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                               family = "curved_gaussian", scalar = curved_gaussian_param,
                               max_iter = max_iter, max_val = max_val,
                               return_path = F, cores = ncores,
                               verbose = F)
save.image("../results/darmanis.RData")

init <- eSVD::initialization(dat_impute, family = "curved_gaussian", k = k, max_val = max_val,
                             scalar = curved_gaussian_param)
curved_gaussian_fit <- eSVD::fit_factorization(dat_impute, u_mat = init$u_mat, v_mat = init$v_mat,
                                                       family = "curved_gaussian", scalar = curved_gaussian_param,
                                                       max_iter = max_iter, max_val = max_val,
                                                       return_path = F, cores = ncores,
                                                       verbose = F)
save.image("../results/darmanis.RData")
