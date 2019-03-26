rm(list=ls())
load("../results/step4_factorization_spca.RData")
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

combn_mat <- combn(4,2)
for(x in 1:ncol(combn_mat)){
  i1 <- combn_mat[1,x]; i2 <- combn_mat[2,x]
  xlim <- range(res_our$u_mat[,i1])
  ylim <- range(res_our$u_mat[,i2])
  order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
  name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
                "Myelin-forming (5)", "Mature (6)")

  png(paste0("../figure/experiment/Week35_marques_latent_spca_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
  par(mfrow = c(2,3), mar = c(4,4,4,0.5))
  for(i in 1:6){
    idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
    plot(res_our$u_mat[-idx,i1], res_our$u_mat[-idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
         main = name_vec[i], xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
           col = rgb(1,1,1), cex = 1.5)

    points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
           col = col_idx[idx], cex = 1.5)
  }
  graphics.off()
}

#############

# let's try a 3d plot
car::scatter3d(x = res_our$u_mat[,1], y = res_our$u_mat[,2], z = res_our$u_mat[,3],
               surface=FALSE, groups = as.factor(cell_type_vec))


################

# plot the raw data
hclust_res <- stats::hclust(stats::dist(t(dat_impute)))

# find out the major cell type partitions
cell_type_vec_raw <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.character(sapply(cell_type_vec_raw, function(x){substr(x,1,2)}))
order_vec <- c(6,5,1,4,2,3)
row_idx <- unlist(lapply(1:length(order_vec), function(i){
  which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
}))
dat2 <- dat_impute[row_idx,hclust_res$order]
cell_type_vec_raw <- cell_type_vec_raw[row_idx]
cell_type_vec <- cell_type_vec[row_idx]
cell_types <- unique(cell_type_vec)

# find out the subtype partitions
for(i in 1:length(order_vec)){
  idx <- which(cell_type_vec == cell_types[i])
  tmp_raw <- cell_type_vec_raw[idx]
  idx2 <- order(tmp_raw)
  cell_type_vec_raw[idx] <- cell_type_vec_raw[idx[idx2]]
  dat2[idx,] <- dat2[idx[idx2], ]
}

# find out where the draw the lines
idx <- which(sapply(1:(length(cell_type_vec)-1), function(x){
  cell_type_vec[x] != cell_type_vec[x+1]
}))
idx_raw <- which(sapply(1:(length(cell_type_vec)-1), function(x){
  cell_type_vec_raw[x] != cell_type_vec_raw[x+1]
}))

stopifnot(all(idx %in% idx_raw))
idx_raw <- idx_raw[!idx_raw %in% idx]
idx <- idx/nrow(dat2)
idx_raw <- idx_raw/nrow(dat2)

png("../figure/experiment/Week35_marques_raw.png", height = 250, width = 1800, res = 300, units = "px")
par(mar = rep(0.5, 4))
.plot_singlecell(t(dat2))

for(i in 1:length(idx)){
  lines(x = rep(idx[i], 2), y = c(0, 1), col = "white", lwd = 2)
  lines(x = rep(idx[i], 2), y = c(0, 1), col = "black", lwd = 2, lty = 2)
}

for(i in 1:length(idx_raw)){
  lines(x = rep(idx_raw[i], 2), y = c(0, 1), col = "white", lwd = 1)
  lines(x = rep(idx_raw[i], 2), y = c(0, 1), col = "black", lwd = 1, lty = 3)
}

graphics.off()

#######################

idx <- which.min(quality_vec)
pred_mat <- 1/(res_list[[idx]]$u_mat %*% t(res_list[[idx]]$v_mat))
mat <- cbind(pred_mat[missing_idx], dat_impute[missing_idx])
pca_res <- stats::princomp(mat)
pca_vec <- pca_res$loadings[,1]

mat_svd <- cbind(pred_naive[missing_idx], dat_impute[missing_idx])
pca_res_svd <- stats::princomp(mat_svd)
pca_vec_svd <- pca_res_svd$loadings[,1]

png("../figure/experiment/Week35_marques_diganostic.png", height = 1000, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2))
plot(NA, asp = T, xlim = range(pred_mat[missing_idx]), ylim = range(dat_impute[missing_idx]),
     xlab = "Predicted value", ylab = "Observed value",
     main = "Prediction of missing values\nvia curved Gaussian (Our method)", cex.main = 0.8)

x_val <- seq(1, 1e5, length.out = 100)
y_val_top <- sapply(x_val, function(x){stats::qnorm(0.9, mean = x, sd = x/scalar_val)})
y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.1, mean = x, sd = x/scalar_val)})

polygon(c(x_val, rev(x_val)), c(y_val_top, rev(y_val_bottom)), col = rgb(1,0,0,0.2),
        border = NA, density = 30, angle = -45)

points(pred_mat[missing_idx], dat_impute[missing_idx], col = rgb(0,0,0,0.2), pch = 16)

lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(x_val, y_val_top, col = "red", lwd = 1, lty = 2)
lines(x_val, y_val_bottom, col = "red", lwd = 1, lty = 2)
lines(c(-1e10,1e10)*pca_vec[1], c(-1e10,1e10)*pca_vec[2], col = "blue", lwd = 2, lty = 2)

plot(NA, asp = T, xlim = range(pred_mat[missing_idx]), ylim = range(dat_impute[missing_idx]),
     xlab = "Predicted value", ylab = "Observed value",
     main = "Prediction of missing values\nvia fixed-value Gaussian (SVD)",
     cex.main = 0.8)

sd_val <- sd(pred_naive[missing_idx]-dat_impute[missing_idx])
x_val <- seq(1, 1e5, length.out = 100)
y_val_top <- sapply(x_val, function(x){stats::qnorm(0.9, mean = x, sd = sd_val)})
y_val_bottom <- sapply(x_val, function(x){stats::qnorm(0.1, mean = x, sd = sd_val)})

polygon(c(x_val, rev(x_val)), c(y_val_top, rev(y_val_bottom)), col = rgb(1,0,0,0.2),
        border = NA, density = 30, angle = -45)

points(pred_naive[missing_idx], dat_impute[missing_idx], pch = 16, col = rgb(0,0,0,0.2))

lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)
lines(x_val, y_val_top, col = "red", lwd = 1, lty = 2)
lines(x_val, y_val_bottom, col = "red", lwd = 1, lty = 2)
lines(c(-1e10,1e10)*pca_vec_svd[1], c(-1e10,1e10)*pca_vec_svd[2], col = "blue", lwd = 2, lty = 2)
graphics.off()

############################

# plot the cell sub-types


# find out the major cell type partitions
cell_type_vec_raw <- as.character(marques$cell.info$cell.type[cell_idx])
combn_mat <- utils::combn(3,2)

alpha_val <- 0.5
col_vec <- c(rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(227/255,73/255,86/255,alpha_val), #red
  rgb(238/255,204/255,17/255,alpha_val), #goldenrod
rgb(129/255,199/255,124/255,alpha_val), #green
rgb(100/255,100/255,200/255,alpha_val), # purple
rgb(40/255,225/255,201/255,alpha_val) #turqouise
)

png("../figure/experiment/Week35_marques_NF.png", height = 1500, width = 2000, res = 300, units = "px")
idx <- grep("NF", cell_type_vec_raw)
len <- length(unique(cell_type_vec_raw[idx]))
par(mfrow = c(len,3))
for(j in 1:len){
  for(i in 1:ncol(combn_mat)){
    i1 <- combn_mat[1,i]; i2 <- combn_mat[2,i]
    idx_j <- which(cell_type_vec_raw[idx] == unique(cell_type_vec_raw[idx])[j])

    if(i==2) main_vec <- paste0("Newly formed sub-type ", j, "\n(", length(idx_j), ")") else main_vec = ""
    plot(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8),
         xlim = range(c(0, res_our$u_mat[idx,i1])), ylim = range(c(0, res_our$u_mat[idx,i2])),
         main = main_vec,
         xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = rgb(1,1,1), cex = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = col_vec[j], cex = 1.5)
  }
}
graphics.off()


png("../figure/experiment/Week35_marques_MF.png", height = 1500, width = 2000, res = 300, units = "px")
idx <- grep("MF", cell_type_vec_raw)
len <- length(unique(cell_type_vec_raw[idx]))
par(mfrow = c(len,3))
for(j in 1:len){
  for(i in 1:ncol(combn_mat)){
    i1 <- combn_mat[1,i]; i2 <- combn_mat[2,i]
    idx_j <- which(cell_type_vec_raw[idx] == unique(cell_type_vec_raw[idx])[j])

    if(i==2) main_vec <- paste0("Myelin-forming sub-type ", j, "\n(", length(idx_j), ")") else main_vec = ""
    plot(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8),
         xlim = range(c(0, res_our$u_mat[idx,i1])), ylim = range(c(0, res_our$u_mat[idx,i2])),
         main = main_vec,
         xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = rgb(1,1,1), cex = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = col_vec[j], cex = 1.5)
  }
}
graphics.off()


png("../figure/experiment/Week35_marques_MO_1.png", height = 2250, width = 2000, res = 300, units = "px")
idx <- grep("MO", cell_type_vec_raw)
len <- length(unique(cell_type_vec_raw[idx]))
par(mfrow = c(3,3))
for(j in c(2,1,6)){
  for(i in 1:ncol(combn_mat)){
    i1 <- combn_mat[1,i]; i2 <- combn_mat[2,i]
    idx_j <- which(cell_type_vec_raw[idx] == unique(cell_type_vec_raw[idx])[j])

    if(i==2) main_vec <- paste0("Mature sub-type ", j, "\n(", length(idx_j), ")") else main_vec = ""
    plot(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8),
         xlim = range(c(0, res_our$u_mat[idx,i1])), ylim = range(c(0, res_our$u_mat[idx,i2])),
         main = main_vec,
         xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = rgb(1,1,1), cex = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = col_vec[j], cex = 1.5)
  }
}
graphics.off()

png("../figure/experiment/Week35_marques_MO_2.png", height = 2250, width = 2000, res = 300, units = "px")
idx <- grep("MO", cell_type_vec_raw)
len <- length(unique(cell_type_vec_raw[idx]))
par(mfrow = c(3,3))
for(j in c(3,4,5)){
  for(i in 1:ncol(combn_mat)){
    i1 <- combn_mat[1,i]; i2 <- combn_mat[2,i]
    idx_j <- which(cell_type_vec_raw[idx] == unique(cell_type_vec_raw[idx])[j])

    if(i==2) main_vec <- paste0("Mature sub-type ", j, "\n(", length(idx_j), ")") else main_vec = ""
    plot(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8),
         xlim = range(c(0, res_our$u_mat[idx,i1])), ylim = range(c(0, res_our$u_mat[idx,i2])),
         main = main_vec,
         xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = rgb(1,1,1), cex = 2)

    points(res_our$u_mat[idx[idx_j],i1], res_our$u_mat[idx[idx_j],i2], pch = 16,
           col = col_vec[j], cex = 1.5)
  }
}
graphics.off()


###########################

# plot the SVD embedding
svd_res <- svd(dat_impute)
u_mat <- svd_res$u[,1:2]%*%diag(svd_res$d[1:2])

alpha_val <- 1
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
plot(u_mat[,1], u_mat[,2], pch = 16, col = col_idx, asp = T) #looks garbage
plot(res_our$u_mat[,1], res_our$u_mat[,2], pch = 16, col = col_idx, asp = T)

set.seed(10)
tsne_res <- tsne::tsne(dat_impute, perplexity = 30, k = 2)

png("../figure/experiment/Week35_marques_tsne.png", height = 1200, width = 1200, res = 300, units = "px")
plot(tsne_res[,1], tsne_res[,2], pch = 16, col = col_idx, asp = T,
     xlab = "Latent dimension 1", ylab = "Latent dimension 2")
graphics.off()
