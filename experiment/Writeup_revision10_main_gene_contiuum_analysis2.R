rm(list=ls())
load("../results/tmp.RData")

cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_factor <- 1

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(zz1[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list,
                               verbose = T, upscale_factor = upscale_factor,
                               reduction_percentage = 0.2,
                               squared = T)
esvd_curves$lineages

####################################

cell_type_vec <- as.character(marques$cell.info$cell.type)
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

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
plotting_order_svd <- c(2,3,1,4)

#####################################3

# fused both sets of trajectory times together

# step 1: for cells in both trajectories -- smartly average the pseudotimes together
mat_list <- .extract_pseudotimes(esvd_curves)
idx <- intersect(mat_list[[1]][,1], mat_list[[2]][,1])
plot(NA, xlim = c(0, 5), ylim = c(0, 5), asp = T)
for(i in idx){
  idx1 <- which(mat_list[[1]][,1] == i)
  idx2 <- which(mat_list[[2]][,1] == i)
  points(mat_list[[1]][idx1, 2], mat_list[[2]][idx2, 2],
         col = col_vec2_svd[cluster_labels[i]],
         pch = 16)
}

# let's just show the "good" cells because 5000 is too many anyways
ideal_cells <- .select_ideal_cells(mat_list)
plot(ideal_cells[,"pseudotime"], col = col_vec2_svd[cluster_labels[ideal_cells[,"cell_idx"]]], pch = 16)

# step 2: determine the relation between all the remaining cells and cells common to both trajectories
trajectory1_mat <- .append_trajectory_specific_cells(pseudotime_mat = mat_list[[1]], cell_mat = ideal_cells,
                                  remaining_idx = which(cluster_labels %in% c(7,6,5)))


other_idx <- which(cluster_labels %in% c(9))
quantile(mat_list[[2]][which(mat_list[[2]][,1] %in% other_idx),2], probs = seq(0,1,length.out = 11))
