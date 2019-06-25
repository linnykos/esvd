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
title(ylab = "(Rescaled) Expression of\n2nd eigen-gene",
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
title(ylab = "(Rescaled) Expression of\n1st eigen-gene",
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
