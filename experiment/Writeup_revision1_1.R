rm(list=ls())
library(xtable)
ncores <- 1

labels_vec <- c("../../raw_data/Lingxue_attachment/Cleaned/Baron_human1/Baron_human1_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_human2/Baron_human2_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_human3/Baron_human3_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_human4/Baron_human4_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse1/Baron_mouse1_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse2/Baron_mouse2_label.csv",
                "../../raw_data/Lingxue_attachment/Cleaned/Darmanis/Darmanis_label.csv")
dat_vec <- c("../../raw_data/Lingxue_attachment/Cleaned/Baron_human1/Baron_human1_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_human2/Baron_human2_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_human3/Baron_human3_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_human4/Baron_human4_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse1/Baron_mouse1_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Baron_mouse2/Baron_mouse2_expr.csv",
             "../../raw_data/Lingxue_attachment/Cleaned/Darmanis/Darmanis_expr.csv")


preprocess_data <- function(dat, label_vec){
  tab <- table(label_vec)
  idx <- which(tab/sum(tab) <= 0.05)

  dat <- dat[which(!label_vec %in% names(tab[idx])),]
  label_vec <- label_vec[which(!label_vec %in% names(tab[idx]))]

  # next remove genes
  tmp <- apply(dat, 2, function(x){sum(x != 0)})
  dat <- dat[,which(tmp >= floor(nrow(dat)/100))]
  dat <- as.matrix(dat)

  # just for now, half the number of cells
  idx <- sample(1:nrow(dat), round(nrow(dat)/2))
  dat <- dat[idx,]
  label_vec <- label_vec[idx]

  # next, follow the esvd preprocessing:

  # try a series of SPCAs
  k <- 5
  lvls <- 5
  v_seq <- exp(seq(log(1), log(log(ncol(dat))/5), length.out = lvls))
  res_list <- vector("list", lvls)

  spca_func <- function(i){
    res <- PMA::SPC(dat, sumabsv = v_seq[i], K = k, trace = F)
    print(paste0("Finished SPC for level ", i))
    res
  }

  print("Starting SPCA")
  res_list <- foreach::"%dopar%"(foreach::foreach(i = 1:lvls), spca_func(i))

  print("Starting DESCEND")
  res_descend <- descend::runDescend(t(dat), n.cores = ncores)

  res_hvg <- descend::findHVG(res_descend, threshold = 50)
  length(res_hvg$HVG.genes)

  print("Compiling results")
  tmp_spca_mat <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
  idx_spca <- which(tmp_spca_mat[,2] == min(tmp_spca_mat[which(tmp_spca_mat[,2] >= 50),2]))
  if(length(idx_spca) == 0) idx_spca <- nrow(tmp_spca_mat)

  idx1 <- sort(unlist(apply(res_list[[idx_spca]]$v, 2, function(x){which(x != 0)})))
  idx2 <- which(colnames(dat) %in% res_hvg$HVG.genes)
  idx <- sort(unique(c(idx1, idx2)))

  list(dat = dat[,idx], label_vec = label_vec, spca_list = res_list, res_hvg = res_hvg)
}


dat_preprocessed <- vector("list", 7)

save.image("../experiment/tmp.RData")

for(i in 1:7){
  print(paste0("Starting data ", i))
  labels <- read.csv(labels_vec[i])
  dat <- read.csv(dat_vec[i])
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  label_vec <- labels$cluster

  dat_preprocessed[[i]] <- preprocess_data(dat, label_vec)

  save.image("../experiment/Writeup_revision1_1.RData")
}

###################

png("../figure/experiment/revision1_data_subsampled.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(2,4), mar = c(0.5, 0.5, 4, 0.5))
for(i in 1:7){
  #rescale each row
  dat_impute <- dat_preprocessed[[i]]$dat
  reweight_factor <- rowSums(dat_impute)
  dat_impute <- t(sapply(1:nrow(dat_impute), function(x){dat_impute[x,]/reweight_factor[x]}))
  # dat_impute <- log(dat_impute+1)
  dat_impute <- dat_impute * 10/mean(dat_impute)

  #simple svd
  svd_res <- svd(dat_impute)
  plot(svd_res$u[,1], svd_res$u[,2], asp = T, col = as.factor(dat_preprocessed[[i]]$label_vec), pch = 16,
       main = i, xaxt = "n", yaxt = "n")
}
graphics.off()

#############

# function to create the neighborhood graph


# i <- 1
# labels <- read.csv(labels_vec[i])
# dat <- read.csv(dat_vec[i])
# dat <- dat[,-1]
# dim(dat)
#
# zz <- apply(dat, 2, function(x){sum(x != 0)})
# length(which(zz >= floor(nrow(dat)/100)))
# length(which(zz >= floor(nrow(dat)/100)))/length(zz)
#
# head(labels)
# zz <- table(labels$cluster)
# tmp <- data.frame(Cluster = names(zz), Number = as.integer(zz))
# xtable::xtable(tmp)
