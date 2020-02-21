rm(list=ls())
load("../results/step4_factorization.RData")

quality_vec

# diagnostic of eSVD
idx <- which.min(quality_vec)
mat_list <- lapply(1:length(res_list[[idx]]), function(j){
  pred_mat <- 1/(res_list[[idx]][[j]]$u_mat %*% t(res_list[[idx]][[j]]$v_mat))

  cbind(dat_impute[missing_idx_list[[j]]], pred_mat[missing_idx_list[[j]]])
})
mat <- do.call(rbind, mat_list)

pca_res <- stats::prcomp(mat, center = F, scale = F)

x_val <- seq(1, 1e5, length.out = 100)
y_val_top <- sapply(x_val, function(x){stats::qnorm(0.9, mean = x, sd = x/scalar_val)})
y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.1, mean = x, sd = x/scalar_val)})

png(filename = "../figure/experiment/Revision_writeup2_main_diagnostic.png",
    height = 1750, width = 1750, res = 300,
    units = "px")
plot(NA, asp = T, pch = 16,
     xlim = c(0, max(mat)),
     ylim = c(0, max(mat)),
     main = paste0("eSVD embedding:\nCurved Gaussian (alpha = ", scalar_val, ")"),
     xlab = "Predicted missing value", ylab = "Observed value",
     cex.lab = 1.25, cex.axis = 1.25)

polygon(c(x_val, rev(x_val)), c(y_val_top, rev(y_val_bottom)), col = rgb(1,0,0,0.2),
        border = NA, density = 30, angle = -45)
points(mat[,2], mat[,1], pch = 16,
       col = rgb(0,0,0,0.2))

lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(c(-1e10,1e10)*pca_res$rotation[2,1], c(-1e10,1e10)*pca_res$rotation[1,1],
      col = "blue", lwd = 2, lty = 2)
lines(x_val, y_val_top, col = "red", lwd = 2, lty = 2)
lines(x_val, y_val_bottom, col = "red", lwd = 2, lty = 2)

rad <- 2/5*max(mat[,1])
ang <- as.numeric(acos(abs(c(1,0) %*% pca_res$rotation[,1])))
radian_seq <- seq(0, ang, length.out = 100)
x_circ <- rad * cos(radian_seq)
y_circ <- rad * sin(radian_seq)
lines(x_circ, y_circ, lty = 2)
text(x = rad, y = 2/5*rad, pos = 4, label = paste0(round(ang * 180/pi, 1), " degrees"))
graphics.off()

#############################################


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


mean_vec <- t(sapply(1:13, function(i){
  idx <- which(cluster_labels == i)
  colMeans(res_our$u_mat[idx,])
}))

png(filename = "../figure/experiment/Revision_writeup2_main_2dplots_1.png",
    height = 1750, width = 1750, res = 300,
    units = "px")
plot(res_our$u_mat[,1], res_our$u_mat[,2], asp = T, pch = 16, col = col_vec3[cluster_labels],
     xlab = "Latent dimension 1", ylab = "Latent dimension 2", main = "eSVD results")
text(mean_vec[,1], mean_vec[,2], labels = as.character(col_info$order), col = "black", cex = 1, font = 2)
graphics.off()

png(filename = "../figure/experiment/Revision_writeup2_main_2dplots_2.png",
    height = 1750, width = 1750, res = 300,
    units = "px")
plot(res_our$u_mat[,1], res_our$u_mat[,3], asp = T, pch = 16, col = col_vec3[cluster_labels],
     xlab = "Latent dimension 1", ylab = "Latent dimension 3", main = "eSVD results")
text(mean_vec[,1], mean_vec[,3], labels = as.character(col_info$order), col = "black", cex = 1, font = 2)
graphics.off()

png(filename = "../figure/experiment/Revision_writeup2_main_2dplots_3.png",
    height = 1750, width = 1750, res = 300,
    units = "px")
plot(res_our$u_mat[,2], res_our$u_mat[,3], asp = T, pch = 16, col = col_vec3[cluster_labels],
     xlab = "Latent dimension 2", ylab = "Latent dimension 2", main = "eSVD results")
text(mean_vec[,2], mean_vec[,3], labels = as.character(col_info$order), col = "black", cex = 1, font = 2)
graphics.off()
