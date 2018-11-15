rm(list=ls())
library(singlecell)
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100,25)/100,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)

.data_generator <- function(cell_pop, gene_pop,
                            distr_func = function(x){stats::rnorm(1, 4/x, sd = 2/x)},
                            n_each = 50, d_each = 120, sigma = 0.05){

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

  res <- .reparameterize(cell_mat, gene_mat)
  cell_mat <- res$X; gene_mat <- res$Y

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- distr_func(max(gram_mat[i,j], 1e-4))
    }
  }

  obs_mat[obs_mat < 0] <- 0

  list(dat = obs_mat, cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = gram_mat, n_each = n_each, d_each = d_each,
       h = h, g = g, k = k, A = res$A)
}

.compute_true_density <- function(cell_pop, gene_pop, grid_size,
                                  sigma = 0.05){
  #discretize the cell line
  cell_mat <- matrix(NA, nrow = 0, ncol = 2)
  for(i in 1:nrow(cell_pop)){
    tmp <- t(sapply(seq(0, 1, length.out = 100), function(x){
      x*cell_pop[i,(1:2)] + (1-x)*cell_pop[i,(3:4)]
    }))
    cell_mat <- rbind(cell_mat, tmp)
  }

  #construct the gene_mat
  g <- nrow(gene_pop)
  gene_mat <- do.call(rbind, lapply(1:g, function(x){
    pos <- stats::runif(200)
    cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(200, sd = sigma),
          pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(200, sd = sigma))
  }))

  pred_mat <- cell_mat %*% t(gene_mat)
  svd_res <- svd(pred_mat)
  cell_mat <- svd_res$u[,1:2] %*% diag(sqrt(svd_res$d[1:2]))
  gene_mat <- svd_res$v[,1:2] %*% diag(sqrt(svd_res$d[1:2]))

  res <- .reparameterize(cell_mat, gene_mat)
  cell_mat <- res$X; gene_mat <- res$Y

  svd_res <- svd(cell_mat); cell_mat <- svd_res$u[,1:2] %*% diag(svd_res$d[1:2])
  spacing <- 0.25
  xrange <- c(floor(min(cell_mat[,1])/spacing)*spacing, ceiling(max(cell_mat[,1])/spacing)*spacing)
  yrange <- c(floor(min(cell_mat[,2])/spacing)*spacing, ceiling(max(cell_mat[,2])/spacing)*spacing)

  xseq <- seq(xrange[1], xrange[2], length.out = grid_size)
  yseq <- seq(yrange[1], yrange[2], length.out = grid_size)
  mat <- matrix(NA, nrow = length(yseq), ncol = length(xseq))

  colnames(mat) <- xseq; rownames(mat) <- yseq

  for(i in 1:length(xseq)){
    if(i %% floor(length(xseq)/10) == 0) cat('*')

    for(j in 1:length(yseq)){
      d <- min(apply(cell_mat, 1, function(x){
        .l2norm(x - c(xseq[i], yseq[j]))
      }))
      mat[j,i] <- dnorm(d, sd = sigma)
    }
  }

  mat <- mat[nrow(mat):1,]

  list(density_mat = mat, cell_mat = cell_mat)
}

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green

#############################

n_seq <- round(seq(10, 50, length.out = 4))
d_seq <- round(seq(20, 100, length.out = 4))
res_list <- vector("list", 4)

for(i in 1:4){
  set.seed(10)
  res <- .data_generator(cell_pop, gene_pop, n_each = n_seq[i], d_each = d_seq[i])
  dat <- res$dat

  # real analysis
  init <- .initialization(dat, family = "gaussian", max_val = 10)
  res_list[[i]] <- .fit_factorization(dat, init$u_mat, init$v_mat,
                            max_val = 5, family = "gaussian", verbose = T,
                            max_iter = 10, reparameterize = T,
                            return_path = F)
}

set.seed(10)
res <- .compute_true_density(cell_pop, gene_pop, 101)
mat <- res$density_mat
mat <- mat[nrow(mat):1, ncol(mat):1]
rownames(mat) <- -as.numeric(rownames(mat))
colnames(mat) <- -as.numeric(colnames(mat))
res$cell_mat <- -res$cell_mat
# image(.rotate(mat))
# plot(res$cell_mat[,1], res$cell_mat[,2])

#identify all the high-probability regions
idx <- which(mat > quantile(mat, probs = 0.95), arr.ind = T)
#translate into coordinates
idx[,1] <- as.numeric(rownames(mat))[idx[,1]]
idx[,2] <- as.numeric(colnames(mat))[idx[,2]]
idx <- idx[,c(2,1)]

#assign each point to clusters
cluster_labels <- as.numeric(apply(idx, 1, function(x){
  i <- which.min(abs(apply(res$cell_mat, 1, function(y){.l2norm(x-y)})))
  floor((i-1)/100)+1
}))
# plot(idx[,1], idx[,2], asp = T, col = col_vec[cluster_labels])

#rescale tmp
r <- range(as.numeric(colnames(mat)))
idx[,1] <- (idx[,1]-r[1])/diff(r)
r <- range(as.numeric(rownames(mat)))
idx[,2] <- (idx[,2]-r[1])/diff(r)

lineages <- list(Lineage1 = c(1,2,3), Lineage2 = c(1,2,4))
b <- quantile(dist(idx[which(cluster_labels == 1),]), probs = 0.05)
pop_curves <- .get_curves(idx, cluster_labels = cluster_labels,
                      lineages, b = 0.005)

png("../figure/simulation/example_trajectories.png", height = 960, width = 2500, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,0.1))

image(.rotate(mat), col = grDevices::heat.colors(100, alpha = 0.5),
      ylab = "X[,2]", xlab = "X[,1]", asp = diff(range(as.numeric(rownames(mat))))/diff(range(as.numeric(colnames(mat)))),
      main = "Population trajectory", cex.lab = 1.25,
      axes = F)
contour(.rotate(mat), add = T, drawlabels = F, col = rgb(0,0,0,0.5), lwd = 1,
        levels = quantile(mat, probs = c(0.25,0.5,0.75)))

xaxs <- as.numeric(colnames(mat)); xaxs <- round(xaxs[round(seq(1,length(xaxs),length.out = 5))],1)
axis(1, at=seq(0,1,length.out = 5), labels = F, las=2)
text(seq(0,1,length.out = 5), par("usr")[3] - 0.1, labels = xaxs, srt = 0, pos = 1,
     xpd = T)
yaxs <- as.numeric(rownames(mat)); yaxs <- round(yaxs[round(seq(1,length(yaxs),length.out = 5))],1)
axis(2, at=seq(0,1,length.out = 5), labels = rev(yaxs), las=2)

for(j in 1:length(pop_curves)){
  lines(pop_curves[[j]], lwd = 2)
}

# plot the centers of each cluster
for(i in 1:4){
  vec <- colMeans(idx[which(cluster_labels == i),])
  points(vec[1], vec[2], pch = 21, cex = 2, bg = col_vec[i])
}

# HOT FIX
for(i in c(1,2,4)){
  tmp <- res_list[[i]]$u_mat; svd_res <- svd(tmp); tmp <- svd_res$u[,1:2] %*% diag(svd_res$d[1:2])

  sig <- ifelse(mean(tmp[1:n_seq[i],2]) > mean(tmp[(3*n_seq[i]+1):(4*n_seq[i]),2]), -1, 1)
  tmp[,2] <- sig*tmp[,2]
  plot(tmp[,1], tmp[,2], xlim = range(as.numeric(colnames(mat))),
       ylim = range(as.numeric(rownames(mat))),
       pch = 16, col = col_vec[rep(1:4, each = n_seq[i])],
       asp = T,
       xlab = "X[,1]", ylab = "X[,2]", axes = F, cex.lab = 1.25,
       main = paste0("Estimated trajectory\n(n = ", 4*n_seq[i], ")"))

  xaxs <- as.numeric(colnames(mat)); xaxs <- round(xaxs[round(seq(1,length(xaxs),length.out = 5))],1)
  axis(1, at = xaxs, labels = F, las=2)
  text(xaxs, par("usr")[3] - 0.1, labels = xaxs, srt = 0, pos = 1,
       xpd = T)
  yaxs <- as.numeric(rownames(mat)); yaxs <- round(yaxs[round(seq(1,length(yaxs),length.out = 5))],1)
  axis(2, at= rev(yaxs), labels = rev(yaxs), las=2)


  lineages <- list(Lineage1 = c(1,2,3), Lineage2 = c(1,2,4))
  b <- quantile(dist(tmp[1:n_seq[i],]), probs = 0.05)
  curves <- .get_curves(tmp, cluster_labels = rep(1:4, each = n_seq[i]),
                        lineages, b = b)
  for(j in 1:length(curves)){
    lines(curves[[j]], lwd = 2)
  }
}
graphics.off()
