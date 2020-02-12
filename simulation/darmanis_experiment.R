rm(list=ls())
library(simulation)
library(eSVD)

##############################

# load the data
#####
# labels_vec <- c("/raid6/Kevin/singlecell_data/Cleaned/Darmanis/Darmanis_label.csv")
# dat_vec <- c("/raid6/Kevin/singlecell_data/Cleaned/Darmanis/Darmanis_expr.csv")
#
# preprocess_data <- function(dat, label_vec){
#   tab <- table(label_vec)
#   idx <- which(tab/sum(tab) <= 0.05)
#
#   dat <- dat[which(!label_vec %in% names(tab[idx])),]
#   label_vec <- label_vec[which(!label_vec %in% names(tab[idx]))]
#
#   # next remove genes
#   tmp <- apply(dat, 2, function(x){sum(x != 0)})
#   dat <- dat[,which(tmp >= floor(nrow(dat)/100))]
#   dat <- as.matrix(dat)
#
#   # just for now, half the number of cells
#   idx <- sample(1:nrow(dat), round(nrow(dat)/2))
#   dat <- dat[idx,]
#   label_vec <- label_vec[idx]
#
#   # next, follow the esvd preprocessing:
#
#   # try a series of SPCAs
#   k <- 5
#   lvls <- 5
#   v_seq <- exp(seq(log(1), log(log(ncol(dat))/5), length.out = lvls))
#   res_list <- vector("list", lvls)
#
#   spca_func <- function(i){
#     res <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
#     print(paste0("Finished SPC for level ", i))
#     res
#   }
#
#   print("Starting SPCA")
#   res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))
#
#   print("Starting DESCEND")
#   res_descend <- descend::runDescend(t(dat), n.cores = ncores)
#
#   res_hvg <- descend::findHVG(res_descend, threshold = 50)
#   length(res_hvg$HVG.genes)
#
#   print("Compiling results")
#   tmp_spca_mat <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
#   idx_spca <- which(tmp_spca_mat[,2] == min(tmp_spca_mat[which(tmp_spca_mat[,2] >= 50),2]))
#   if(length(idx_spca) == 0) idx_spca <- nrow(tmp_spca_mat)
#
#   idx1 <- sort(unlist(apply(res_list[[idx_spca]]$v, 2, function(x){which(x != 0)})))
#   idx2 <- which(colnames(dat) %in% res_hvg$HVG.genes)
#   idx <- sort(unique(c(idx1, idx2)))
#
#   dat <- dat[,idx]
#
#   #rescale each row
#   reweight_factor <- rowSums(dat)
#   dat <- t(sapply(1:nrow(dat), function(x){dat[x,]/reweight_factor[x]}))
#   dat <- dat * 10/mean(dat)
#
#   list(dat = dat, label_vec = label_vec, spca_list = res_list, res_hvg = res_hvg)
# }
# labels <- read.csv(labels_vec[1])
# dat <- read.csv(dat_vec[1])
# dat <- dat[,-1]
# dat <- as.matrix(dat)
# label_vec <- labels$cluster
#
# dat_preprocessed <- preprocess_data(dat, label_vec)
#############

load("../experiment/Writeup_revision1_1.RData")
dat <- dat_preprocessed[[7]]$dat
reweight_factor <- rowSums(dat)
dat <- t(sapply(1:nrow(dat), function(x){dat[x,]/reweight_factor[x]}))
dat <- dat * 10/mean(dat)
var_all <- ls()
var_all <- var_all[!var_all %in% "dat"]
rm(list=var_all)

###

paramMat <- cbind(c(1:3, rep(4,5)), c(rep(1,3), 1:5), 3, 50)
colnames(paramMat) <- c("distr", "param", "k", "max_iter")
trials <- 5
ncores <- 20
doMC::registerDoMC(cores = ncores)

save.image("../experiment/darmanis_experiment_data.RData")

#################################

rule <- function(vec){
  dat
}

criterion <- function(dat, vec, y){
  set.seed(10*y)

  n <- nrow(dat); d <- ncol(dat)
  missing_idx <- rbind(do.call(rbind, (lapply(1:n, function(x){
    cbind(x, sample(1:d, 2))
  }))), do.call(rbind, (lapply(1:d, function(x){
    cbind(sample(1:n, 2), d)
  }))))

  dat_NA <- dat
  for(i in 1:nrow(missing_idx)){
    dat_NA[missing_idx[i,1], missing_idx[i,2]] <- NA
  }
  missing_idx <- which(is.na(dat_NA))
  missing_val <- dat[missing_idx]

  set.seed(10)
  if(vec["distr"] == 1){
    init <- eSVD::initialization(dat, family = "gaussian", k = vec["k"], max_val = 100)
    tmp <- eSVD::fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "gaussian",
                                   max_iter = vec["max_iter"], max_val = 100,
                                   return_path = F, cores = ncores,
                                   verbose = F)


    pred_mat <- tmp$u_mat %*% t(tmp$v_mat)
    pred_val <- pred_mat[missing_idx]

  } else if(vec["distr"] == 2){
    init <- eSVD::initialization(dat, family = "exponential", k = vec["k"], max_val = -100)
    tmp <- eSVD::fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "exponential",
                                   max_iter = vec["max_iter"], max_val = -100,
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- tmp$u_mat %*% t(tmp$v_mat)
    pred_val <- -1/pred_mat[missing_idx]


  } else if(vec["distr"] == 3){
    init <- eSVD::initialization(dat, family = "poisson", k = vec["k"], max_val = 100)
    tmp <- eSVD::fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "poisson",
                                   max_iter = vec["max_iter"], max_val = 100,
                                   return_path = F, cores = ncores,
                                   verbose = F)

    pred_mat <- tmp$u_mat %*% t(tmp$v_mat)
    pred_val <- exp(pred_mat[missing_idx])

  } else {
    init <- eSVD::initialization(dat, family = "curved_gaussian", k = vec["k"], max_val = 100,
                                 scalar = vec["param"])
    tmp <- eSVD::fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                   family = "curved_gaussian", scalar = vec["param"],
                                   max_iter = vec["max_iter"], max_val = 100,
                                   return_path = F, cores = ncores,
                                   verbose = T)

    pred_mat <- tmp$u_mat %*% t(tmp$v_mat)
    pred_val <- 1/pred_mat[missing_idx]
  }

  list(fit = tmp, pred_val = pred_val, missing_val = missing_val)
}

##

##########################################


res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../experiment/darmanis_experiment_tmp.RData",
                                        verbose = T)

save.image("../experiment/darmanis_experiment.RData")
