# load(paste0("../results/step5_clustering", suffix, ".RData"))
load(paste0("../results/old_results/step5_clustering_spca.RData"))

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha)) #orange
}
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)

######################################################

# diagnostic of eSVD
idx <- which.min(quality_vec)
pred_mat <- 1/(res_list[[idx]]$u_mat %*% t(res_list[[idx]]$v_mat))
tmp <- cbind(-pred_mat[missing_idx], dat_impute[missing_idx])
pca_res <- princomp(tmp)

x_val <- seq(1, 1e5, length.out = 100)
y_val_top <- sapply(x_val, function(x){stats::qnorm(0.75, mean = x, sd = x/scalar_val)})
y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.25, mean = x, sd = x/scalar_val)})

png("../figure/main/esvd_diagnostic.png", height = 1200, width = 1100, res = 300, units = "px")
plot(NA, asp = T, pch = 16,
     xlim = c(0, max(c(-pred_mat[missing_idx], dat_impute[missing_idx]))),
     ylim = c(0, max(c(-pred_mat[missing_idx], dat_impute[missing_idx]))),
     main = "eSVD embedding:\nMatrix-completion diagnostic",
     xlab = "Predicted value", ylab = "Observed but masked value",
     cex.lab = 1.25, cex.axis = 1.25)

polygon(c(x_val, rev(x_val)), c(y_val_top, rev(y_val_bottom)), col = rgb(1,0,0,0.2),
        border = NA, density = 30, angle = -45)
points(-pred_mat[missing_idx], dat_impute[missing_idx], pch = 16,
       col = rgb(0,0,0,0.2))

lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
      col = "blue", lwd = 2, lty = 2)
lines(x_val, y_val_top, col = "red", lwd = 2, lty = 2)
lines(x_val, y_val_bottom, col = "red", lwd = 2, lty = 2)

rad <- 500
ang <- as.numeric(acos(c(1,0) %*% pca_res$loadings[,1]))
radian_seq <- seq(0, ang, length.out = 100)
x_circ <- rad * cos(radian_seq)
y_circ <- rad * sin(radian_seq)
lines(x_circ, y_circ, lty = 2)
text(x = 500, y = 200, pos = 4, label = paste0(round(ang * 180/pi, 1), " degrees"))
graphics.off()

######################################

#diagnostic of SVD
tmp <- cbind(pred_naive[missing_idx], dat_impute[missing_idx])
pca_res <- princomp(tmp)

x_val <- seq(1, 1e5, length.out = 100)
sd_val <- stats::sd(pred_naive[-missing_idx] - dat_impute[-missing_idx])
y_val_top <- sapply(x_val, function(x){stats::qnorm(0.9, mean = x, sd = sd_val)})
y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.1, mean = x, sd = sd_val)})

png("../figure/main/svd_diagnostic.png", height = 1200, width = 1100, res = 300, units = "px")
plot(NA, asp = T, pch = 16,
     xlim = c(0, max(c(pred_naive[missing_idx], dat_impute[missing_idx]))),
     ylim = c(0, max(c(pred_naive[missing_idx], dat_impute[missing_idx]))),
     main = "SVD embedding:\nMatrix-completion diagnostic",
     xlab = "Predicted value", ylab = "Observed but masked value",
     cex.lab = 1.25, cex.axis = 1.25)

polygon(c(x_val, rev(x_val)), c(y_val_top, rev(y_val_bottom)), col = rgb(1,0,0,0.2),
        border = NA, density = 30, angle = -45)
points(pred_naive[missing_idx], dat_impute[missing_idx], pch = 16,
       col = rgb(0,0,0,0.2))

lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
lines(c(-1e10,1e10)*pca_res$loadings[1,1], c(-1e10,1e10)*pca_res$loadings[1,2],
      col = "blue", lwd = 2, lty = 2)
lines(x_val, y_val_top, col = "red", lwd = 2, lty = 2)
lines(x_val, y_val_bottom, col = "red", lwd = 2, lty = 2)

rad <- 500
ang <- as.numeric(acos(c(1,0) %*% pca_res$loadings[,1]))
radian_seq <- seq(0, ang, length.out = 100)
x_circ <- rad * cos(radian_seq)
y_circ <- rad * sin(radian_seq)
lines(x_circ, y_circ, lty = 2)
text(x = 500, y = 200, pos = 4, label = paste0(round(ang * 180/pi, 1), " degrees"))
graphics.off()

##################

svd_res <- svd(dat_impute)
svd_u <- svd_res$u[,1:p] %*% diag(sqrt(svd_res$d[1:p]))

col_vec <- color_func(1)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]
col_name <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

png("../figure/main/svd_preview.png", height = 1200, width = 2200, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(5,6,4,2))

scaling_vec <- -svd_res$v[,2]
scaling_vec <- (scaling_vec-min(scaling_vec))/(max(scaling_vec)-min(scaling_vec))
scaling_vec <- scaling_vec/sum(scaling_vec)
val <- dat_impute %*% scaling_vec
tmp_df <- data.frame(val = val, type = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))})[cluster_labels])

vioplot::vioplot(tmp_df$val[tmp_df$type == 1],
                 tmp_df$val[tmp_df$type == 2],
                 tmp_df$val[tmp_df$type == 3],
                 tmp_df$val[tmp_df$type == 4],
                 tmp_df$val[tmp_df$type == 5],
                 tmp_df$val[tmp_df$type == 6],
                 col = sapply(1:6, function(x){unique(col_info$col_code[which(col_info$level == x)])}),
                 pchMed = 21,
                 colMed = "black", colMed2 = "white",
                 xlab = "", names = rep("", 6))
title(ylab = "Weighted expression based\non 2nd singular vector",
      main = "Average gene expression\nper cell type (SVD)", cex.lab = 1.25)
text(1:6, par("usr")[3]-2,
     srt = -45, xpd = TRUE,
     labels = c("Pdgfra+", "OPC", "COP", "NFOL", "MFOL", "MOL"), cex=1)

par(mar = c(5,4,4,2))
set.seed(10)
idx <- sample(1:nrow(svd_u))
plot(-svd_u[idx,2], svd_u[idx,3], asp = T, pch = 16, col = col_vec[cluster_labels][idx],
     cex = 0.75, xlab = "Latent dimension 2", ylab = "Latent dimension 3",
     main = "SVD embedding\n(Constant-variance Gaussian)",
     cex.lab = 1.25, cex.axis = 1.25)
graphics.off()

####################################

our_curves$lineages
custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))

col_vec <- color_func(1)[c(5, rep(3,2), 3, rep(1,3), rep(4,2), rep(2,2),  rep(5,2))]
col_name <- c("orange", rep("bluish green", 3), rep("yellow", 3), rep("blue", 2), rep("skyblue", 2), rep("orange", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(custom_cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

nat_mat <- res_our$u_mat %*% t(res_our$v_mat)
svd_res <- svd(nat_mat)
scaling_vec <- -svd_res$v[,1]
scaling_vec <- (scaling_vec-min(scaling_vec))/(max(scaling_vec)-min(scaling_vec))
scaling_vec <- scaling_vec/sum(scaling_vec)
val <- dat_impute %*% scaling_vec

# custom tmp_df
type_vec <- sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))})[cluster_labels]
tmp_df <- data.frame(val = val, type = type_vec)
tmp_df <- rbind(tmp_df, data.frame(val = val[cluster_labels == 4], type = 7)) #duplicate MOL shared in both lineages

png("../figure/main/esvd_preview.png", height = 1200, width = 2200, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(5,6,4,2))

vioplot::vioplot(tmp_df$val[tmp_df$type == 1],
                 tmp_df$val[tmp_df$type == 2],
                 tmp_df$val[tmp_df$type == 3],
                 tmp_df$val[tmp_df$type == 4],
                 tmp_df$val[tmp_df$type == 5],
                 tmp_df$val[tmp_df$type == 6],
                 tmp_df$val[tmp_df$type == 7],
                 col = sapply(1:7, function(x){
                   tmp <- table(col_info$col_code[which(col_info$level == x)])
                   names(tmp)[which.max(tmp)]
                 }),
                 pchMed = 21,
                 colMed = "black", colMed2 = "white",
                 xlab = "", names = rep("", 7))
title(ylab = "Weighted expression based\non 1st singular vector",
      main = "Average gene expression\nper cell type (eSVD)", cex.lab = 1.25)
text(1:7, par("usr")[3]-3,
     srt = -45, xpd = TRUE,
     labels = c("Pdgfra+", "OPC", "COP", "NFOL", "MFOL", "MOL (1)", "MOL (2)"), cex=.75)

par(mar = c(5,4,4,2))
set.seed(10)
idx <- sample(1:nrow(res_our$u_mat))
plot(res_our$u_mat[idx,1], res_our$u_mat[idx,3], asp = T, pch = 16, col = col_vec[cluster_labels][idx],
     cex = 0.75, xlab = "Latent dimension 1", ylab = "Latent dimension 3",
     main = "eSVD embedding\n(Curved Gaussian)",
     cex.lab = 1.25, cex.axis = 1.25)
graphics.off()

########################################

# 2D plots

cluster_center <- .compute_cluster_center(res_our$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))
custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))

num_order_vec <- c(1, rep(3,3), rep(5,3), rep(4,2), rep(2,2),  rep(1,2))
col_vec <- color_func(1)[num_order_vec]
col_vec2 <- color_func(0.1)[num_order_vec]
col_vec3 <- color_func(0.3)[num_order_vec]
col_name <- c("yellow", rep("bluish green", 3), rep("orange", 3), rep("blue", 2), rep("skyblue", 2), rep("yellow", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(custom_cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

png(paste0("../figure/main/eSVD_trajectory_2d.png"),
    height = 1000, width = 2700, res = 300, units = "px")
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
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 7)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = col_vec_short[ll], lwd = 5)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
          lty = 3, lwd = 2)
  }
}
graphics.off()

# 3D plots

png(paste0("../figure/main/eSVD_theta45_phi225.png"),
    height = 2000, width = 2000, res = 300, units = "px")
par(mar = c(0,0,4,0))
slingshot_3dplot(res_our$u_mat[,1:3], cluster_labels,
                 bg_col_vec = col_vec2,
                 breaks = seq(0.5, 13.5, by = 1),
                 cluster_center = cluster_center,
                 center_col_vec = col_vec,
                 center_labels = 1:13,
                 curves = NA,
                 pch = 16, main = "eSVD embedding and trajectories\n(Curved Gaussian)",
                 xlab = "", ylab = "",
                 zlab = "Latent dimension 3",
                 theta = 45, phi = 225,
                 xlim = c(-2.45, 0), ylim = c(-1.5,1), zlim = c(-1,0.75))

curves <- our_curves$curves
col_vec_short <- color_func(0.9)[c(5,4)]
for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                  y = curves[[i]]$s[ord, 2],
                  z = curves[[i]]$s[ord, 3],
                  add = T, colkey = F, col = "black", lwd = 2)
  plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                  y = curves[[i]]$s[ord, 2],
                  z = curves[[i]]$s[ord, 3],
                  add = T, colkey = F, col = col_vec_short[i], lwd = 2)
}
graphics.off()

# 3D plots with tubes

our_tube_list <- lapply(1:length(our_curves$curves), function(x){
  s_mat <- our_curves$curves[[x]]$s[our_curves$curves[[x]]$ord,]
  construct_3d_tube(s_mat, radius = our_sd_val$sd_val)
})

png(paste0("../figure/main/eSVD_theta45_phi225_tube.png"),
    height = 2000, width = 2000, res = 300, units = "px")
par(mar = c(0,0,4,0))
slingshot_3dplot(res_our$u_mat[,1:3], cluster_labels,
                 bg_col_vec = col_vec2,
                 breaks = seq(0.5, 13.5, by = 1),
                 cluster_center = cluster_center,
                 center_col_vec = col_vec,
                 center_labels = 1:13,
                 curves = NA,
                 pch = 16, main = "eSVD embedding with uncertainty tubes\n(Curved Gaussian)",
                 xlab = "", ylab = "",
                 zlab = "Latent dimension 3",
                 theta = 45, phi = 225,
                 xlim = c(-2.45, 0), ylim = c(-1.5,1), zlim = c(-1,0.75))

curves <- our_curves$curves
col_vec_short <- color_func(0.9)[c(5,4)]
for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                  y = curves[[i]]$s[ord, 2],
                  z = curves[[i]]$s[ord, 3],
                  add = T, colkey = F, col = "black", lwd = 2)
  plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                  y = curves[[i]]$s[ord, 2],
                  z = curves[[i]]$s[ord, 3],
                  add = T, colkey = F, col = col_vec_short[i], lwd = 2)

  col_mat <- our_tube_list[[i]]$z_mat

  plot3D::surf3D(our_tube_list[[i]]$x_mat,
                 our_tube_list[[i]]$y_mat,
                 our_tube_list[[i]]$z_mat, add = T,
                 colvar = col_mat,
                 col = colorRampPalette(c("white", col_vec_short[i]))(100),
                 breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                 colkey = F)
}
graphics.off()

###########################


svd_res <- svd(dat_impute)
svd_u <- svd_res$u[,1:p] %*% diag(sqrt(svd_res$d[1:p]))

cluster_center <- .compute_cluster_center(svd_u[,1:3], .construct_cluster_matrix(cluster_labels))

num_order_vec <- c(1, rep(3,2), rep(5,6), rep(2,2),  rep(1,2))
col_vec <- color_func(1)[num_order_vec]
col_vec3 <- color_func(0.3)[num_order_vec]
col_name <- c("yellow", rep("bluish green", 2), rep("orange", 6), rep("skyblue", 2), rep("yellow", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

# SVD 2d plots

png(paste0("../figure/main/SVD_trajectory_2d.png"),
    height = 1000, width = 2700, res = 300, units = "px")
par(mfrow = c(1,3))
combn_mat <- combn(3,2)
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]
  plot(x = svd_u[,i], y = svd_u[,j],
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       col = col_vec3[cluster_labels], pch = 16,
       main = ifelse(k == 2, "SVD embedding and trajectories\n(Constant-variance Gaussian)","")
  )

  for(ll in 1:nrow(cluster_center)){
    points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 3, col = "black")
    points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 2, col = col_vec[[ll]])
  }


  curves <- naive_curves$curves
  for(ll in 1:length(curves)){
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 4)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 1.5)
  }
}
graphics.off()


# svd 3D plots

png(paste0("../figure/main/SVD_theta135_phi180.png"),
    height = 2000, width = 2000, res = 300, units = "px")
par(mar = c(1,1,4,1))
slingshot_3dplot(svd_u[,1:3], cluster_labels,
                 bg_col_vec = col_vec2,
                 breaks = seq(0.5, 13.5, by = 1),
                 cluster_center = cluster_center,
                 center_col_vec = col_vec,
                 center_labels = 1:13,
                 curves = NA,
                 pch = 16, main = "SVD embedding and trajectories\n(Constant-variance Gaussian)",
                 xlab = "", ylab = "Latent dimension 2",
                 zlab = "Latent dimension 3",
                 theta = 135, phi = 180)

curves <- naive_curves$curves
col_vec_short <- color_func(0.9)[c(5,4)]
for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                  y = curves[[i]]$s[ord, 2],
                  z = curves[[i]]$s[ord, 3],
                  add = T, colkey = F, col = "black", lwd = 2)
}
graphics.off()

# svd 3D plots with tubes

naive_tube_list <- lapply(1:length(naive_curves$curves), function(x){
  s_mat <- naive_curves$curves[[x]]$s[naive_curves$curves[[x]]$ord,]
  construct_3d_tube(s_mat, radius = naive_sd_val$sd_val)
})
col_vec_short <- color_func(0.9)[c(4)]

png(paste0("../figure/main/SVD_theta135_phi180_tube.png"),
    height = 2000, width = 2000, res = 300, units = "px")
par(mar = c(1,1,4,1))
slingshot_3dplot(svd_u[,1:3], cluster_labels,
                 bg_col_vec = col_vec2,
                 breaks = seq(0.5, 13.5, by = 1),
                 cluster_center = cluster_center,
                 center_col_vec = col_vec,
                 center_labels = 1:13,
                 curves = NA,
                 pch = 16, main = "SVD embedding and trajectories\n(Constant-variance Gaussian)",
                 xlab = "", ylab = "Latent dimension 2",
                 zlab = "Latent dimension 3",
                 theta = 135, phi = 180)

curves <- naive_curves$curves
col_vec_short <- color_func(0.9)[c(4)]
for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                  y = curves[[i]]$s[ord, 2],
                  z = curves[[i]]$s[ord, 3],
                  add = T, colkey = F, col = "black", lwd = 2)
}

i <- 3
col_mat <- naive_tube_list[[i]]$z_mat

plot3D::surf3D(naive_tube_list[[i]]$x_mat,
               naive_tube_list[[i]]$y_mat,
               naive_tube_list[[i]]$z_mat, add = T,
               colvar = col_mat,
               col = colorRampPalette(c("white", col_vec_short[1]))(100),
               breaks = seq(min(col_mat), max(col_mat), length.out = 101),
               colkey = F)
graphics.off()

