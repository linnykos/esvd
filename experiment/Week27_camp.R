rm(list=ls())
load("../../SOUP/data/camp.rda")

dat <- camp$counts

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 10*log(dat + 1)

dim(dat)
dat <- dat[order(camp$cell.info[,2]),]
camp$cell.info <- camp$cell.info[order(camp$cell.info[,2]),]

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

png("../figure/experiment/27_camp_need_imputed.png", height = 800, width = 1500, res = 300, units = "px")
par(mar = rep(0.5,4))
zero_mat2 <- zero_mat
zero_mat2[zero_mat2 == 0] <- 1
zero_mat2[is.na(zero_mat2)] <- 0
image(.rotate(zero_mat2), breaks = c(-0.5,0.5,1.5), col = c("brown1", "gray88"),
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
.plot_singlecell(dat_impute)
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
k <- 2
u_mat <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

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
graphics.off()
