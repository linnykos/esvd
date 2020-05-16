rm(list=ls())
load("../results/step5_trajectory_original.RData")

# count how many zeros
length(which(dat_impute == 0))/prod(dim(dat_impute)) # 36% zeros

######

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}

num_order_vec_svd <- c(5, rep(3,2), c(1,1,1,1,1,1), rep(2,2),  rep(5,2))
col_vec_svd <- color_func(1)[num_order_vec_svd]
col_vec2_svd <- color_func(0.5)[num_order_vec_svd]
col_name_svd <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
order_vec_svd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_svd <- data.frame(name = levels(cell_type_vec),
                           idx = sort(unique(cluster_labels)),
                           order = order_vec_svd,
                           col_name = col_name_svd,
                           col_code = col_vec_svd)
col_info_svd$factor_idx <- as.numeric(as.factor(col_info_svd$col_name))
col_info_svd[,c(5,6)] <- col_info_svd[,c(6,5)]
colnames(col_info_svd)[c(5,6)] <- colnames(col_info_svd)[c(6,5)]
col_info_svd

######
for(i in 1:4){
  cluster_val <- col_info_svd$idx[which(col_info_svd$factor_idx == i)]
  dat_subset <- dat_impute[which(cluster_labels %in% cluster_val),]
  print(length(which(dat_subset == 0))/prod(dim(dat_subset)))
}
# huh... pretty weird. the youngest cell type has the most zeros, but we see in the fits it often has the highest predicted values...
# pretty odd...

for(i in 1:4){
  cluster_val <- col_info_svd$idx[which(col_info_svd$factor_idx == i)]
  dat_subset <- dat_impute[which(cluster_labels %in% cluster_val),]
  print(quantile(dat_subset))
}

for(i in 1:4){
  cluster_val <- col_info_svd$idx[which(col_info_svd$factor_idx == i)]
  dat_subset <- dat_impute[which(cluster_labels %in% cluster_val),]
  print(sd(dat_subset))
}

for(i in 1:4){
  cluster_val <- col_info_svd$idx[which(col_info_svd$factor_idx == i)]
  dat_subset <- dat_impute[which(cluster_labels %in% cluster_val),]
  print(sd(dat_subset[dat_subset != 0]))
}

