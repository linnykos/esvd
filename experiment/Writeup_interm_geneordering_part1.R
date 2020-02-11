rm(list=ls())
load("../results/old_results/step4_factorization_spca.RData")

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_vec <- rep(NA, length(unique(cluster_labels)))
size_vec <- sapply(cluster_group_list, function(x){length(which(cluster_labels %in% x))})
for(i in 1:length(cluster_group_list)){
  upscale_vec[cluster_group_list[[i]]] <- (max(size_vec)/size_vec[i])^(1/2)
}


p <- 3
set.seed(10)
our_curves <- slingshot(res_our$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                              cluster_group_list = cluster_group_list,
                              verbose = T, upscale_vec = upscale_vec)

save.image("../experiment/Writeup_interm_geneordering.RData")
#####

# make sure the curves look reasonable

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(210/255, 198/255, 36/255, alpha)) #darker yellow
}
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)


cluster_center <- .compute_cluster_center(res_our$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))
custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))

num_order_vec <- c(5, rep(3,2), 3, rep(1,3), rep(4,2), rep(2,2),  rep(5,2))
col_vec <- color_func(1)[num_order_vec]
col_vec2 <- color_func(0.1)[num_order_vec]
col_vec3 <- color_func(0.3)[num_order_vec]
col_name <- c("orange", rep("bluish green", 3), rep("yellow", 3), rep("blue", 2), rep("skyblue", 2), rep("orange", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(custom_cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

col_vec_short <- color_func(0.9)[c(6,4)]

par(mfrow = c(1,3))
combn_mat <- combn(3,2)
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]
  plot(x = res_our$u_mat[,i], y = res_our$u_mat[,j],
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       col = col_vec3[cluster_labels], pch = 16,
       main = ifelse(k == 2, "eSVD embedding and trajectories\n(Curved Gaussian)","")
  )

  for(ll in 1:nrow(cluster_center)){
    points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 3, col = "black")
    points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 2, col = col_vec[[ll]])
  }


  curves <- our_curves$curves
  for(ll in 1:length(curves)){
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = col_vec_short[ll], lwd = 5)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
          lty = 3, lwd = 2)
  }
}
