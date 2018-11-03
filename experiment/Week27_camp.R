rm(list=ls())
load("../../SOUPR/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 10*log(dat + 1)

dim(dat)
dat <- dat[order(camp$cell.info[,2]),]
camp$cell.info <- camp$cell.info[order(camp$cell.info[,2]),]

############################

# plot some histograms
dropout_res <- lapply(1:ncol(dat), function(i){
  .em_mixture(dat[,i])
})



# highest mean
# lowest mean

png("../figure/experiment/27_camp_histograms.png", height = 2000, width = 2500, res = 300, units = "px")
par(mfrow = c(2,2), mar = c(4,3,3,1))
zz <- which.max(apply(dat, 1, function(x){length(which(x == 0))})) # most zeros
.hist_augment(dat[zz,], main = paste0("Cell ", zz), xlab = "Value")
zz <- which.min(apply(dat, 1, function(x){length(which(x == 0))})) # least zeros
.hist_augment(dat[zz,], main = paste0("Cell ", zz), xlab = "Value")

zz <- which.min(sapply(dropout_res, function(x){x$class2["mean"]}))
.hist_augment(dat[,zz], main = paste0("Gene ", zz), xlab = "Value")
zz <- which.max(sapply(dropout_res, function(x){x$class2["mean"]}))
.hist_augment(dat[,zz], main = paste0("Gene ", zz), xlab = "Value")
graphics.off()

################

cell_type <- camp$cell.info[,2]
cell_type_coarse <- as.numeric(as.factor(sapply(cell_type, function(x){substr(x, 1, 1)})))
cell_type_numeric <- as.numeric(as.factor(cell_type))

################

res_svd <- svd(dat)

k <- 4
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

row_idx <- unlist(lapply(1:max(cell_type_coarse), function(x){
  which(cell_type_coarse == x)
}))
row_idx_lines <- sapply(1:(max(cell_type_coarse)-1), function(x){
  1-length(which(cell_type_coarse <= x))/length(cell_type_coarse)
})
col_idx <- unlist(lapply(1:max(v_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))
col_idx_lines <- sapply(1:(max(u_clust$cluster)-1), function(x){
  1-length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

dat <- dat[, col_idx]

png("../figure/experiment/27_camp_data.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
.plot_singlecell(dat)
for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}

for(i in col_idx_lines){
  lines(rep(i, 2), c(0,1), lwd = 3, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 3, lty = 2)
}
graphics.off()


################

dropout_mat <- singlecell:::.dropout(dat)
zero_mat <- singlecell:::.find_true_zeros(dropout_mat, num_neighbors = 50)
length(which(zero_mat == 0))/prod(dim(zero_mat))

png("../figure/experiment/27_camp_true_zero.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
zero_mat2 <- zero_mat
zero_mat2[is.na(zero_mat2)] <- 1
image(.rotate(zero_mat2), breaks = c(-0.5,0.5,1.5), col = c("blue3", "gray88"),
      asp = nrow(zero_mat2)/ncol(zero_mat2),
      axes = F, main = "")
for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}

for(i in col_idx_lines){
  lines(rep(i, 2), c(0,1), lwd = 3, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 3, lty = 2)
}
graphics.off()

png("../figure/experiment/27_camp_zero_mat.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
zero_mat2 <- zero_mat
zero_mat2[zero_mat2 == 0] <- 1
zero_mat2[is.na(zero_mat2)] <- 0
image(.rotate(zero_mat2), breaks = c(-0.5,0.5,1.5), col = c("goldenrod1", "gray88"),
      asp = nrow(zero_mat2)/ncol(zero_mat2),
      axes = F, main = "")
for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}

for(i in col_idx_lines){
  lines(rep(i, 2), c(0,1), lwd = 3, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 3, lty = 2)
}
graphics.off()

set.seed(10)
dat_impute <- singlecell:::.scImpute(dat, which(is.na(zero_mat)), Kcluster = 7,
                                     verbose = T, weight = 0.25)

png("../figure/experiment/27_camp_imputed.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
.plot_singlecell(dat_impute, luminosity = F)
for(i in row_idx_lines){
  lines(c(0,1), rep(i, 2), lwd = 3, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 3, lty = 2)
}

for(i in col_idx_lines){
  lines(rep(i, 2), c(0,1), lwd = 3, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 3, lty = 2)
}
graphics.off()

svd_res <- svd(dat_impute)
k <- 4
u_mat <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
# res <- .identification(u_mat, v_mat)
# u_mat <- res$X; v_mat <- res$Y

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green
col_vec2 <- c("coral", "brown", "gray70", "black",
              "cyan", "dodgerblue", "blue")

png("../figure/experiment/27_camp_latent.png", height = 1200, width = 1200, res = 300, units = "px")
plot(u_mat[,1], u_mat[,2],
     xlim = range(c(u_mat[,1], 0)), ylim = range(c(u_mat[,2], 0)),
     col = col_vec[cell_type_coarse], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

png("../figure/experiment/27_camp_latent_detailed.png", height = 1200, width = 1200, res = 300, units = "px")
plot(u_mat[,1], u_mat[,2],
     xlim = range(c(u_mat[,1], 0)), ylim = range(c(u_mat[,2], 0)),
     col = col_vec2[cell_type_numeric], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

cluster_labels <- as.numeric(as.factor(cell_type))
lineages <- list(Lineage1 = c("7", "6", "5", "4", "3", "2", "1"))
cluster_mat <- sapply(1:7, function(x){
  vec <- rep(0, length(cluster_labels))
  vec[which(cluster_labels == x)] <- 1
  vec
})
colnames(cluster_mat) <- 1:7
# res <- .get_lineages(dat, cluster_labels) #LOL doesn't work rip
# lineages <- res$lineages
# lineages[[1]] <- c("1","3","2")
# lineages[[2]] <- c("1","3","4")
# cluster_mat <- res$cluster_mat
curves <- .get_curves_tmp(u_mat, cluster_mat, lineages, reassign = F, b = 0.03)

png("../figure/experiment/27_camp_latent_withcurves.png", height = 1200, width = 1200, res = 300, units = "px")
plot(u_mat[,1], u_mat[,2],
     xlim = range(c(u_mat[,1], 0)), ylim = range(c(u_mat[,2], 0)),
     col = col_vec[cell_type_coarse], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
ord <- curves[[1]]$ord
lines(curves[[1]]$s[ord,1], curves[[1]]$s[ord,2], lwd = 2)

#find the centers
center_points <- sapply(1:max(cluster_labels), function(x){
  idx <- which(cluster_labels == x)
  colMeans(u_mat[idx,])[1:2]
})

#find the corresponding places on the curve
curve_points <- sapply(1:ncol(center_points), function(x){
  idx <- which.min(apply(curves[[1]]$s, 1, function(y){
    .l2norm(y[1:2] - center_points[,x])
  }))
  curves[[1]]$s[idx,]
})

for(i in 1:7){
  points(center_points[1,i], center_points[2,i], pch = 16)
  points(curve_points[1,i], curve_points[2,i], pch = 16)
  lines(c(center_points[1,i], curve_points[1,i]), c(center_points[2,i], curve_points[2,i]),
        lwd = 2, lty = 2)
}

# label_tag <- c(4,2,4,3,4,4,4)
# text(x = center_points[1,], y = center_points[2,], labels = unique(cell_type),
#      pos = rev(label_tag), cex = 0.75)

lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()


######################################
load("../experiment/Week27_camp.RData")
col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green
cell_type <- camp$cell.info[,2]
cell_type_coarse <- as.numeric(as.factor(sapply(cell_type, function(x){substr(x, 1, 1)})))
cell_type_numeric <- as.numeric(as.factor(cell_type))

plot(res_withimpute$u_mat[,1], res_withimpute$u_mat[,2],
     xlim = range(c(res_withimpute$u_mat[,1], 0)), ylim = range(c(res_withimpute$u_mat[,2], 0)),
     col = col_vec[cell_type_coarse], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")

col_vec2 <- c("coral", "brown", "gray70", "black",
              "cyan", "dodgerblue", "blue")

plot(res_withimpute$u_mat[,1], res_withimpute$u_mat[,2],
     xlim = range(c(res_withimpute$u_mat[,1], 0)), ylim = range(c(res_withimpute$u_mat[,2], 0)),
     col = col_vec2[cell_type_numeric], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")

plot(res_withimpute$v_mat[,1], res_withimpute$v_mat[,2],
     xlim = range(c(res_withimpute$v_mat[,1], 0)), ylim = range(c(res_withimpute$v_mat[,2], 0)),
     asp = T, pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors")

#######################
# actually plot the results now

png("../figure/experiment/27_camp_estimates.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(res_withimpute$u_mat[,1], res_withimpute$u_mat[,2],
     xlim = range(c(res_withimpute$u_mat[,1], 0)), ylim = range(c(res_withimpute$u_mat[,2], 0)),
     col = col_vec[cell_type_coarse], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(res_withimpute$v_mat[,1], res_withimpute$v_mat[,2],
     xlim = range(c(res_withimpute$v_mat[,1], 0)), ylim = range(c(res_withimpute$v_mat[,2], 0)),
     asp = T, pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors",
     col = rgb(0,0,0,0.5))
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()


#########

cluster_labels <- as.numeric(as.factor(cell_type))
lineages <- list(Lineage1 = c("7", "6", "5", "4", "2", "1", "3"))
cluster_mat <- sapply(1:7, function(x){
  vec <- rep(0, length(cluster_labels))
  vec[which(cluster_labels == x)] <- 1
  vec
})
colnames(cluster_mat) <- 1:7
# res <- .get_lineages(dat, cluster_labels) #LOL doesn't work rip
# lineages <- res$lineages
# lineages[[1]] <- c("1","3","2")
# lineages[[2]] <- c("1","3","4")
# cluster_mat <- res$cluster_mat
curves <- .get_curves_tmp(res_withimpute$u_mat, cluster_mat, lineages, reassign = F, b = 5)

#find the centers
center_points <- sapply(1:max(cluster_labels), function(x){
  idx <- which(cluster_labels == x)
  colMeans(res_withimpute$u_mat[idx,])[1:2]
})

#find the corresponding places on the curve
curve_points <- sapply(1:ncol(center_points), function(x){
  idx <- which.min(apply(curves[[1]]$s, 1, function(y){
    .l2norm(y[1:2] - center_points[,x])
  }))
  curves[[1]]$s[idx,]
})

png("../figure/experiment/27_camp_estimates_withcurves.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(res_withimpute$u_mat[,1], res_withimpute$u_mat[,2],
     xlim = range(c(res_withimpute$u_mat[,1], 0)), ylim = range(c(res_withimpute$u_mat[,2], 0)),
     col = col_vec[cell_type_coarse], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
ord <- curves[[1]]$ord
lines(curves[[1]]$s[ord,1], curves[[1]]$s[ord,2], lwd = 2)

lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

for(i in 1:7){
  points(center_points[1,i], center_points[2,i], pch = 16)
  points(curve_points[1,i], curve_points[2,i], pch = 16)
  lines(c(center_points[1,i], curve_points[1,i]), c(center_points[2,i], curve_points[2,i]),
        lwd = 2, lty = 2)
}

plot(res_withimpute$v_mat[,1], res_withimpute$v_mat[,2],
     xlim = range(c(res_withimpute$v_mat[,1], 0)), ylim = range(c(res_withimpute$v_mat[,2], 0)),
     asp = T, pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors",
     col = rgb(0,0,0,0.5))
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()


png("../figure/experiment/27_camp_estimates_preandpost.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(res_withimpute$u_mat[,1], res_withimpute$u_mat[,2],
     xlim = range(c(res_withimpute$u_mat[,1], 0)), ylim = range(c(res_withimpute$u_mat[,2], 0)),
     col = col_vec[cell_type_coarse], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
ord <- curves[[1]]$ord
lines(curves[[1]]$s[ord,1], curves[[1]]$s[ord,2], lwd = 2)

lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

for(i in 1:7){
  points(center_points[1,i], center_points[2,i], pch = 16)
  points(curve_points[1,i], curve_points[2,i], pch = 16)
  lines(c(center_points[1,i], curve_points[1,i]), c(center_points[2,i], curve_points[2,i]),
        lwd = 2, lty = 2)
}

cluster_labels <- as.numeric(as.factor(cell_type))
lineages <- list(Lineage1 = c("7", "6", "5", "4", "3", "2", "1"))
cluster_mat <- sapply(1:7, function(x){
  vec <- rep(0, length(cluster_labels))
  vec[which(cluster_labels == x)] <- 1
  vec
})
colnames(cluster_mat) <- 1:7
curves <- .get_curves_tmp(u_mat, cluster_mat, lineages, reassign = F, b = 0.03)

plot(u_mat[,1], u_mat[,2],
     xlim = range(c(u_mat[,1], 0)), ylim = range(c(u_mat[,2], 0)),
     col = col_vec[cell_type_coarse], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
ord <- curves[[1]]$ord
lines(curves[[1]]$s[ord,1], curves[[1]]$s[ord,2], lwd = 2)

#find the centers
center_points <- sapply(1:max(cluster_labels), function(x){
  idx <- which(cluster_labels == x)
  colMeans(u_mat[idx,])[1:2]
})

#find the corresponding places on the curve
curve_points <- sapply(1:ncol(center_points), function(x){
  idx <- which.min(apply(curves[[1]]$s, 1, function(y){
    .l2norm(y[1:2] - center_points[,x])
  }))
  curves[[1]]$s[idx,]
})

for(i in 1:7){
  points(center_points[1,i], center_points[2,i], pch = 16)
  points(curve_points[1,i], curve_points[2,i], pch = 16)
  lines(c(center_points[1,i], curve_points[1,i]), c(center_points[2,i], curve_points[2,i]),
        lwd = 2, lty = 2)
}

# label_tag <- c(4,2,4,3,4,4,4)
# text(x = center_points[1,], y = center_points[2,], labels = unique(cell_type),
#      pos = rev(label_tag), cex = 0.75)

lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

