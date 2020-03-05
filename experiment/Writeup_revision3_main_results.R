rm(list=ls())
load("../results/step5_clustering.RData")
scalar_val

# first the diagnostic plots
# first up, fixed variance gaussian
all_mat <- do.call(rbind, naive_res_list)
plot(all_mat[,1], all_mat[,2], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5), asp = T)

# next, for curved Gaussian
all_mat <- do.call(rbind, res_list)
plot(all_mat[,1], all_mat[,2], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5), asp = T)
lines(c(-1e6,1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

pca_res <- stats::prcomp(all_mat[,2:1], center = F, scale = F)

x_val <- seq(1, 1e5, length.out = 100)
y_val_top <- sapply(x_val, function(x){stats::qnorm(0.9, mean = x, sd = x/scalar_val)})
y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.1, mean = x, sd = x/scalar_val)})

lines(c(-1e10,1e10)*pca_res$rotation[2,1], c(-1e10,1e10)*pca_res$rotation[1,1],
      col = "blue", lwd = 2, lty = 2)

#########################################

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
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})


col_vec <- color_func(1)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]
col_name <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
order_vec <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                         order = order_vec,
                       col_name = col_name,
                       col_code = col_vec)
num_order_vec <- c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))
col_vec3 <- color_func(0.3)[num_order_vec]

our_curves$lineages

mean_vec <- t(sapply(1:13, function(i){
  idx <- which(cluster_labels == i)
  colMeans(res_our$u_mat[idx,])
}))

combn_mat <- combn(3,2)
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]
  plot(x = res_our$u_mat[,i], y = res_our$u_mat[,j],
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       col = col_vec3[cluster_labels], pch = 16,
       main = ifelse(k == 2, "eSVD embedding and trajectories\n(Curved Gaussian)","")
  )

  curves <- our_curves$curves
  for(ll in 1:length(curves)){
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
          lty = 3, lwd = 2)
  }

  text(mean_vec[,i], mean_vec[,j], labels = as.character(col_info$order), col = "black", cex = 1, font = 2)
}


