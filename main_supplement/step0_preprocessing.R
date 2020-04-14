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

  save.image(paste0("../results/baron_step0_screening", suffix, ".RData"))
}

rm(list=c("dat", "labels", "label_vec", "tab", "idx", "tmp", "k", "lvls", "v_seq", "res_list", "spca_mat",
          "target_var", "spca_idx", "res_descend", "res_hvg", "descend_idx", "gene_idx", "dat_impute", "reweight_factor"))

print(paste0(Sys.time(), ": Finished preprocessing"))
save.image(paste0("../results/baron_step0_screening", suffix, ".RData"))
