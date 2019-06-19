load(paste0("../results/step5_clustering", suffix, ".RData"))

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

###################

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha)) #orange
}

svd_res <- svd(dat_impute)
svd_u <- svd_res$u[,1:p] %*% diag(sqrt(svd_res$d[1:p]))

col_vec <- color_func(1)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]
col_name <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
cbind(levels(cell_type_vec), sort(unique(cluster_labels)),
      sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))}),
      col_name, col_vec)

png("../figure/main/svd_embedding_23.png", height = 1200, width = 1100, res = 300, units = "px")
plot(svd_u[,2], svd_u[,3], asp = T, pch = 16, col = col_vec[cluster_labels],
     cex = 0.75, xlab = "Latent dimension 2", ylab = "Latent dimension 3",
     main = "SVD embedding\n(Constant-variance Gaussian)",
     cex.lab = 1.25, cex.axis = 1.25)
graphics.off()

##################

plot(res_our$u_mat[,1], res_our$u_mat[,2], asp = T, pch = 16, col = col_vec[cluster_labels],
     cex = 0.75, xlab = "Latent dimension 2", ylab = "Latent dimension 3",
     main = "eSVD embedding\n(Curved Gaussian)",
     cex.lab = 1.25, cex.axis = 1.25)

#########
#under construction


dat <- res_our$u_mat[,1:d]

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
alpha_val <- 0.2
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]

num_cell <- length(unique(cell_type_vec))

combn_mat <- combn(3,2)
for(x in 1:ncol(combn_mat)){
  i1 <- combn_mat[1,x]; i2 <- combn_mat[2,x]
  xlim <- range(dat[,i1])
  ylim <- range(dat[,i2])
  order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
  name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
                "Myelin-forming (5)", "Mature (6)")

  png(paste0("../figure/main/marques_latent_our_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
  par(mfrow = c(2,3), mar = c(4,4,4,0.5))
  for(i in 1:6){
    idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
    plot(dat[-idx,i1], dat[-idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
         main = name_vec[i], xlab = paste0("Our latent dimension ", i1),
         ylab = paste0("Our latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = rgb(1,1,1), cex = 1.5)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = col_idx[idx], cex = 1.5)
  }
  graphics.off()

}
