rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat <- zeisel$counts
dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 1000*log(dat + 1)
#
sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})

sum_vec <- colSums(dat)

quantile_sum <- quantile(sum_vec, prob = 0.95)

idx <- which(sum_vec <= quantile_sum)
dat <- dat[,idx]

sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})

sum_vec <- colSums(dat)

idx <- which(colnames(dat) %in% zeisel$select.genes)
n <- nrow(dat)

png("../figure/experiment/sparsity_3.png", height = 1200, width = 1200, res = 300, units = "px")
col_vec <- rep(rgb(0,0,0,0.1), n)
col_vec[idx] <- rgb(0.803, 0.156, 0.211, 0.5)
plot(sparsity_vec, sum_vec, ylim = c(0, quantile_sum), col = col_vec, pch = 16,
     xlab = "Sparsity (Number of non-zeroes)", ylab = "Summation of counts")
graphics.off()

###############

dat <- dat[,idx]
res_svd <- svd(dat)

png("../figure/experiment/scree_3.png", height = 1200, width = 1600, res = 300, units = "px")
plot(res_svd$d, xlim = c(1, 50), pch = 16,
     xlab = "Rank", ylab = "Singular value")
graphics.off()

u_mat <- res_svd$u[,1:8] %*% diag(res_svd$d[1:8])
v_mat <- res_svd$v[,1:8] %*% diag(res_svd$d[1:8])

png("../figure/experiment/latent_3.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

plot(NA, xlim = range(c(u_mat[,1], 0)),
     ylim = range(c(u_mat[,2], 0)), xlab = "U[,1]", ylab = "U[,2]", main = "Latent vectors for cells")
lines(c(-1e10, 1e10), rep(0, 2), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(rep(0, 2), c(-1e10, 1e10), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
points(u_mat[,1], u_mat[,2], pch = 16, col = rgb(0, 0, 0, 0.1))

plot(NA, xlim = range(c(v_mat[,1], 0)),
     ylim = range(c(v_mat[,2], 0)), xlab = "V[,1]", ylab = "V[,2]", main = "Latent vectors for genes")
lines(c(-1e10, 1e10), rep(0, 2), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(rep(0, 2), c(-1e10, 1e10), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
points(v_mat[,1], v_mat[,2], pch = 16, col = rgb(0, 0, 0, 0.1))

graphics.off()

#####################

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

k_vec <- 2:10
u_clust_vec <- sapply(k_vec, function(x){
  set.seed(10)
  kmeans(u_mat_spherical, centers = x, iter.max = 100, nstart = 10)$tot.withinss
})
v_clust_vec <- sapply(k_vec, function(x){
  set.seed(10)
  kmeans(v_mat_spherical, centers = x, iter.max = 100, nstart = 10)$tot.withinss
})

png("../figure/experiment/kmeans_3.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

plot(k_vec, u_clust_vec,  pch = 16, xlab = "K", ylab = "Total within ss", main = "Clustering for cells")
plot(k_vec, v_clust_vec,pch = 16, xlab = "K", ylab = "Total within ss", main = "Clustering for genes")

graphics.off()

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 6, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 6, iter.max = 100, nstart = 10)

col_clean <- c(rgb(0,0,0,0.5), rgb(0.803, 0.156, 0.211, 0.5), rgb(0.415, 0.564, 0.792, 0.5),
               rgb(0.537, 0.521, 0.901, 0.5), rgb(0.584, 0.858, 0.564, 0.5), rgb(0.988, 0.607, 0.113, 0.5))

png("../figure/experiment/latent_clust_3.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

plot(NA, xlim = range(c(u_mat[,1], 0)),
     ylim = range(c(u_mat[,2], 0)), xlab = "U[,1]", ylab = "U[,2]", main = "Latent vectors for cells")
lines(c(-1e10, 1e10), rep(0, 2), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(rep(0, 2), c(-1e10, 1e10), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
points(u_mat[,1], u_mat[,2], pch = 16, col = col_clean[u_clust$cluster])


plot(NA, xlim = range(c(v_mat[,1], 0)),
     ylim = range(c(v_mat[,2], 0)), xlab = "V[,1]", ylab = "V[,2]", main = "Latent vectors for genes")
lines(c(-1e10, 1e10), rep(0, 2), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(rep(0, 2), c(-1e10, 1e10), col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
points(v_mat[,1], v_mat[,2], pch = 16, col = col_clean[v_clust$cluster])

graphics.off()

##############

u_tab <- table(zeisel$cell.info[,2], u_clust$cluster)
u_tab2 <- t(apply(u_tab, 1, function(x){x/sum(x)}))

uv_tab <- sapply(1:6, function(x){
  sapply(1:6, function(y){
    idx1 <- which(u_clust$cluster == x)
    idx2 <- which(v_clust$cluster == y)
    mean(u_mat[idx1,] %*% t(v_mat[idx2,]))
  })
})

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

png("../figure/experiment/table_3.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(1,4,4,0))
col_ramp <- colorRampPalette(c("white", rgb(0.803, 0.156, 0.211)))(100)

y_tic <- seq(-.5, 6.5, by = 1)*1/6
x_tic <- seq(-.5, 5.5, by = 1)*1/5
image(.rotate(u_tab2), col = col_ramp, asp = 7/5, axes = F,
      xlim = c(min(x_tic)-.05, max(x_tic)+.05),
      ylim = c(min(y_tic)-.05, max(y_tic)+.05), main = "Cluster table for Cells\n(Count)")
for(y in y_tic){
  lines(range(x_tic), rep(y, 2))
}
for(x in x_tic){
  lines(rep(x, 2), range(y_tic))
}
text(par("usr")[3] - 0.1, seq(0,1,length.out = 7), adj = 1, labels = rev(rownames(u_tab)), xpd = TRUE,
     cex = 0.9)
text(seq(0, 1, length.out = 6), par("usr")[1] + 0.3, labels = paste0("C", 1:6), xpd = TRUE,
     cex = 0.9)

x_vec <- seq(0,1,length.out = 6)
y_vec <- seq(0,1,length.out = 7)
for(i in 1:length(x_vec)){
  for(j in 1:length(y_vec)){
    text(x_vec[i], 1-y_vec[j], label = u_tab[j,i])
  }
}

y_tic <- seq(-.5, 5.5, by = 1)*1/5
x_tic <- seq(-.5, 5.5, by = 1)*1/5
image(.rotate(uv_tab), col = col_ramp, asp = T, axes = F,
      xlim = c(min(x_tic)-.05, max(x_tic)+.05),
      ylim = c(min(y_tic)-.05, max(y_tic)+.05), main = "Cluster table for Cells x Genes\n(Mean inner product)")
for(y in y_tic){
  lines(range(x_tic), rep(y, 2))
}
for(x in x_tic){
  lines(rep(x, 2), range(y_tic))
}
text(par("usr")[3], seq(0,1,length.out = 6), adj = 1, labels = rev(paste0("G", 1:6)), xpd = TRUE,
     cex = 0.9)
text(seq(0, 1, length.out = 6), par("usr")[1] + 0.03, labels = paste0("C", 1:6), xpd = TRUE,
     cex = 0.9)

x_vec <- seq(0,1,length.out = 6)
y_vec <- seq(0,1,length.out = 6)
for(i in 1:length(x_vec)){
  for(j in 1:length(y_vec)){
    text(x_vec[i], 1-y_vec[j], label = round(uv_tab[j,i], 1))
  }
}

graphics.off()

#######################

rank <- 8
dat_approx <- res_svd$u[,1:rank] %*% diag(res_svd$d[1:rank]) %*% t(res_svd$v[,1:rank])

.slice_matrix <- function(dat, num.bins){
  x.cuts <- seq(from = min(dat[,1]), to = max(dat[,1]), length = num.bins + 1)
  y.cuts <- seq(from = min(dat[,2]), to = max(dat[,2]), length = num.bins + 1)

  index.x <- cut(dat[,1], x.cuts, include.lowest = TRUE)
  index.y <- cut(dat[,2], y.cuts, include.lowest = TRUE)

  mat <- tapply(dat[,1], list(index.x, index.y), base::length)
  mat[is.na(mat)] <- 0

  mat
}

plot_mat <- cbind(as.numeric(dat_approx), as.numeric(dat))
plot_mat <- plot_mat[intersect(intersect(which(plot_mat[,1] <= 1), which(plot_mat[,2] <= 1)),
                               intersect(which(plot_mat[,1] >= 0), which(plot_mat[,2] >= 0))),]
count_mat <- .slice_matrix(plot_mat, 20)
non_zero_mat <- as.numeric(count_mat)
non_zero_mat <- non_zero_mat[non_zero_mat != 0]

png("../figure/experiment/residual_3.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

plot(as.numeric(dat_approx), as.numeric(dat), col = rgb(0,0,0,0.1), pch = 16, asp = T,
     xlab = "Predicted value", ylab = "Observed value", xlim = range(as.numeric(dat_approx)),
     ylim = range(as.numeric(dat)))
mean_vec <- seq(0, 50, length.out = 1e3)
sd_top <- sapply(mean_vec, function(x){qexp(p = 0.9, rate = 1/x)})
sd_bot <- sapply(mean_vec, function(x){qexp(p = 0.1, rate = 1/x)})
lines(mean_vec, mean_vec, col = rgb(0.803, 0.156, 0.211), lwd = 2)
lines(mean_vec, sd_top, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(mean_vec, sd_bot, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)

plot(plot_mat[,1], plot_mat[,2], col = rgb(0,0,0,0.005), pch = 16, asp = T,
     xlab = "Predicted value", ylab = "Observed value", xlim = range(plot_mat[,1]),
     ylim = range(plot_mat[,2]))
lines(mean_vec, mean_vec, col = rgb(0.803, 0.156, 0.211), lwd = 2)
lines(mean_vec, sd_top, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(mean_vec, sd_bot, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)

contour(seq(from = min(plot_mat[,1]), to = max(plot_mat[,1]), length = nrow(count_mat)),
        seq(from = min(plot_mat[,2]), to = max(plot_mat[,2]), length = ncol(count_mat)),
        count_mat, levels = quantile(non_zero_mat, probs = rev(seq(0, 1, length.out = 50)))[c(2:5,10,15,20,25)], add = T, axes = F, drawlabels = F,
        col = rgb(0.584, 0.858, 0.564), lwd = 2)
graphics.off()

#################

# simulation

set.seed(10)
tol <- 1e-5
mean_mat <- pmax(dat_approx, tol)
obs_mat <- mean_mat
for(i in 1:nrow(mean_mat)){
  for(j in 1:ncol(mean_mat)){
    obs_mat[i,j] <- rexp(1, rate = 1/mean_mat[i,j])
  }
}

obs_svd <- svd(obs_mat)
pred_mean_mat <- obs_svd$u[,1:8] %*% diag(obs_svd$d[1:8]) %*% t(obs_svd$v[,1:8])

plot_mat <- cbind(as.numeric(pred_mean_mat), as.numeric(obs_mat))
plot_mat <- plot_mat[intersect(intersect(which(plot_mat[,1] <= 1), which(plot_mat[,2] <= 1)),
                               intersect(which(plot_mat[,1] >= 0), which(plot_mat[,2] >= 0))),]
count_mat <- .slice_matrix(plot_mat, 20)
non_zero_mat <- as.numeric(count_mat)
non_zero_mat <- non_zero_mat[non_zero_mat != 0]

png("../figure/experiment/simulation_3.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))

plot(as.numeric(pred_mean_mat), as.numeric(obs_mat), col = rgb(0,0,0,0.1), pch = 16, asp = T,
     xlab = "Predicted value", ylab = "Observed value", xlim = range(as.numeric(pred_mean_mat)),
     ylim = range(as.numeric(obs_mat)))
mean_vec <- seq(0, 50, length.out = 1e3)
sd_top <- sapply(mean_vec, function(x){qexp(p = 0.9, rate = 1/x)})
sd_bot <- sapply(mean_vec, function(x){qexp(p = 0.1, rate = 1/x)})
lines(mean_vec, mean_vec, col = rgb(0.803, 0.156, 0.211), lwd = 2)
lines(mean_vec, sd_top, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(mean_vec, sd_bot, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)

plot(plot_mat[,1], plot_mat[,2], col = rgb(0,0,0,0.005), pch = 16, asp = T,
     xlab = "Predicted value", ylab = "Observed value", xlim = range(plot_mat[,1]),
     ylim = range(plot_mat[,2]))
lines(mean_vec, mean_vec, col = rgb(0.803, 0.156, 0.211), lwd = 2)
lines(mean_vec, sd_top, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)
lines(mean_vec, sd_bot, col = rgb(0.803, 0.156, 0.211), lwd = 2, lty = 2)

contour(seq(from = min(plot_mat[,1]), to = max(plot_mat[,1]), length = nrow(count_mat)),
        seq(from = min(plot_mat[,2]), to = max(plot_mat[,2]), length = ncol(count_mat)),
        count_mat, levels = quantile(non_zero_mat, probs = rev(seq(0, 1, length.out = 50)))[c(2:5,10,15,20,25)], add = T, axes = F, drawlabels = F,
        col = rgb(0.584, 0.858, 0.564), lwd = 2)
graphics.off()

