rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- cbind(50, 120, 0.01, 150, 3, 2, 10)
colnames(paramMat) <- c("n", "d", "sigma", "total", "k", "scalar", "max_val")
trials <- 50

################

# setup
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)

.data_generator <- function(cell_pop, gene_pop,
                            n_each = 50, d_each = 100, sigma = 0.05,
                            scalar = 2, total = 150){
  #construct the cell information
  h <- nrow(cell_pop)
  cell_mat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = sigma),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = sigma))
  }))
  n <- nrow(cell_mat)
  k <- ncol(cell_mat)

  # construct the gene information
  g <- nrow(gene_pop)
  gene_mat <- do.call(rbind, lapply(1:g, function(x){
    pos <- stats::runif(d_each)
    cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(d_each, sd = sigma),
          pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(d_each, sd = sigma))
  }))
  d <- nrow(gene_mat)

  # form observations
  gram_mat <- cell_mat %*% t(gene_mat) #natural parameter
  svd_res <- svd(gram_mat)
  cell_mat <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
  gene_mat <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

  res <- singlecell:::.reparameterize(cell_mat, gene_mat)
  cell_mat <- res$u_mat; gene_mat <- res$v_mat

  extra_weight <- rep(1, nrow(cell_mat))
  pred_mat <- 1/(cell_mat %*% t(gene_mat))
  pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
    pred_mat[x,] * extra_weight[x]
  }))

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- rnorm(1, pred_mat[i,j], pred_mat[i,j]/scalar)
    }
  }

  obs_mat[obs_mat < 0] <- 0
  obs_mat2 <- round(exp(obs_mat*14)-1)
  # length(which(obs_mat2 > 5000))/prod(dim(obs_mat2))
  obs_mat2[obs_mat2 > 5000] <- 5000
  # quantile(obs_mat2)

  # now do something more dramatic with dropout
  obs_mat3 <- obs_mat2
  .dropped_indices <- function(x, total){
    vec <- 1:length(x)
    samp <- sample(vec, size = total, replace = T, prob = x)
    setdiff(vec, unique(samp))
  }

  total_vec <- rep(total, nrow(obs_mat3))
  for(i in 1:nrow(obs_mat3)){
    idx <- .dropped_indices(obs_mat[i,], total = total_vec[i])
    obs_mat3[i,idx] <- 0
  }

  list(dat = obs_mat3, dat_nodropout = obs_mat2,
       extra_weight = extra_weight,
       cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = cell_mat %*% t(gene_mat), n_each = n_each, d_each = d_each,
       h = h, g = g, k = k)
}

set.seed(1)
vec <- paramMat[1,]
obj <- .data_generator(cell_pop, gene_pop, n_each = vec["n"], d_each = vec["d"],
                       sigma = vec["sigma"], scalar = vec["scalar"], total = vec["total"])

# add columns to the data to "screen out"
set.seed(1)
dat <- obj$dat
d_add <- 100
dat_new <- matrix(0, nrow = nrow(dat), ncol = d_add)
mean_val <- mean(dat[dat!=0]); sd_val <- sd(dat[dat!=0])
for(i in 1:d_add){
  dat_new[,i] <- rnorm(nrow(dat_new), mean_val, sd_val)
}
dat_new[dat_new < 0] <- 0
dat_new[sample(prod(dim(dat_new)), round(0.7*prod(dim(dat_new))))] <- 0
dat_new[,sample(ncol(dat_new), ceiling(0.8*ncol(dat_new)))] <- 0
dat_new <- cbind(dat, dat_new)
ord <- sample(ncol(dat_new))
dat_new2 <- dat_new[,ord]

png("../figure/illustration/dat_full.png", height = 1100, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
.plot_singlecell(dat_new2)
graphics.off()

png("../figure/illustration/dat.png", height = 1100, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
dat <- obj$dat
.plot_singlecell(dat)
graphics.off()

####################

load("../results/factorization_results.RData")

png("../figure/illustration/dat_impute.png", height = 1100, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
dat <- res[[1]][[1]]$dat_impute
.plot_singlecell(dat)
graphics.off()

png("../figure/illustration/embedding.png", height = 1200, width = 1000, res = 300, units = "px")
dat <- res[[1]][[1]]$res_our$u_mat
plot(dat[,1], dat[,2], xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                 asp = T, col = rgb(0.7,0.7,0.7), pch = 16)
graphics.off()

png("../figure/illustration/clustering.png", height = 1200, width = 1000, res = 300, units = "px")
dat <- res[[1]][[1]]$res_our$u_mat
plot(dat[,1], dat[,2], xlab = "Latent dimension 1", ylab = "Latent dimension 2",
     asp = T, col = rgb(0.7,0.7,0.7), pch = 16)
graphics.off()


col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green
png("../figure/illustration/clustering.png", height = 1200, width = 1000, res = 300, units = "px")
dat <- res[[1]][[1]]$res_our$u_mat
plot(dat[,1], dat[,2], xlab = "Latent dimension 1", ylab = "Latent dimension 2",
     asp = T, col = col_vec[rep(1:4, each = vec["n"])], pch = 16)
graphics.off()

png("../figure/illustration/lineage.png", height = 1200, width = 1000, res = 300, units = "px")
dat <- res[[1]][[1]]$res_our$u_mat
cluster_labels <- c(1:4)[rep(1:4, each = vec["n"])]
b_est <- .b_estimate(res[[1]][[1]]$cell_mat[,1:2], cluster_labels)
fixed_lineage <- singlecell::slingshot(res[[1]][[1]]$cell_mat[,1:2], cluster_labels, 1, knn = NA,
                                       b = b_est, remove_outlier = F)$lineages
b_est <- .b_estimate(dat[,1:2], cluster_labels)
our_lineage <- .get_curves(dat[,1:2], cluster_labels, fixed_lineage, shrink = 1,
                           thresh =  0.001, max_iter = 15, b = b_est)

plot(dat[,1], dat[,2], xlab = "Latent dimension 1", ylab = "Latent dimension 2",
     asp = T, col = rgb(0.7,0.7,0.7), pch = 16)

for(i in 1:length(our_lineage)){
  ord <- our_lineage[[i]]$ord
  lines(our_lineage[[i]]$s[ord, 1], our_lineage[[i]]$s[ord, 2], lwd = 2,
        col = "black")
}

cluster_mat <- .construct_cluster_matrix(cluster_labels)
centers <- .compute_cluster_center(dat, cluster_mat)
points(centers[,1], centers[,2], col = "white", pch = 16, cex = 1.75)
points(centers[,1], centers[,2], col = col_vec, pch = 16, cex = 1.5)

graphics.off()
