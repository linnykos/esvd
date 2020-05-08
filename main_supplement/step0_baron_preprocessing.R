labels_file_vec <- c("/raid6/Kevin/singlecell_data/Cleaned/Baron_human1/Baron_human1_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_human2/Baron_human2_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_human3/Baron_human3_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_human4/Baron_human4_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse1/Baron_mouse1_label.csv",
                     "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse2/Baron_mouse2_label.csv")
dat_file_vec <- c("/raid6/Kevin/singlecell_data/Cleaned/Baron_human1/Baron_human1_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_human2/Baron_human2_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_human3/Baron_human3_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_human4/Baron_human4_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse1/Baron_mouse1_expr.csv",
                  "/raid6/Kevin/singlecell_data/Cleaned/Baron_mouse2/Baron_mouse2_expr.csv")

preprocessing_list <- vector("list", length(dat_file_vec))
for(i in 1:length(labels_file_vec)){
  print(paste0("Starting dataset ", i))

  labels <- read.csv(labels_file_vec[i])
  dat <- read.csv(dat_file_vec[i])
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  label_vec <- labels$cluster

  ##########

  tab <- table(label_vec)
  idx <- which(tab/sum(tab) <= 0.05)

  dat <- dat[which(!label_vec %in% names(tab[idx])),]
  dat_count <- dat
  label_vec <- label_vec[which(!label_vec %in% names(tab[idx]))]
  label_vec <- as.factor(as.character(label_vec))

  # next remove genes
  zz <- apply(dat, 2, function(x){length(which(x!=0))})
  idx <- which(zz > nrow(dat)/100)
  dat <- dat[,idx]
  dat_count <- dat_count[,idx]

  # normalize
  dat <- t(sapply(1:nrow(dat), function(i){10^4 * dat[i,]/sum(dat[i,])}))
  rownames(dat) <- 1:nrow(dat)

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

  # run sPCA
  spca_mat <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
  idx <- min(which(spca_mat[,2] == ncol(dat)))
  target_var <- spca_mat[idx,3]
  idx <- min(intersect(which(spca_mat[,2] >= 30), which(spca_mat[,3] >= 0.9*target_var)))
  spca_idx <- sort(unlist(apply(res_list[[idx]]$v, 2, function(x){which(x != 0)})))
  spca_hvg <- colnames(dat)[spca_idx]

  # run DESCEND
  res_descend <- descend::runDescend(t(dat_count), n.cores = ncores)
  descend_hvg <- descend::findHVG(res_descend, threshold = 50)$HVG.genes

  idx <- which(colnames(dat) %in% c(descend_hvg, spca_hvg))
  dat <- dat[,idx]

  preprocessing_list[[i]] <- list(dat_impute = dat, label_vec = label_vec,
                                  vst_hvg = vst_hvg, spca_hvg = spca_hvg,
                                  spca_mat = spca_mat)

  save.image(paste0("../results/step0_baron_preprocessing", suffix, ".RData"))
}

print(paste0(Sys.time(), ": Finished preprocessing"))
save.image(paste0("../results/step0_baron_screening", suffix, ".RData"))
print(warnings())
