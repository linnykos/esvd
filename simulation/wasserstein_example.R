rm(list=ls())
library(singlecell)

.data_generator <- function(distr_func = function(x){stats::rnorm(1, 4/x, sd = 2/x)},
                            n_each = 50, d_each = 120, sigma = 0.05,
                            multiplier = 1){

  #construct the cell information
  cell_pop <- multiplier*matrix(c(4,10, 25,100,
                                  60,80, 25,100,
                                  40,10, 60,80,
                                  60,80, 100,25)/10,
                                nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  cell_mat_org <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = sigma),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = sigma))
  }))
  n <- nrow(cell_mat_org)
  k <- ncol(cell_mat_org)

  # construct the gene information
  gene_pop <- multiplier*matrix(c(20, 90, 25, 100,
                          90,20, 100,25)/10, nrow = 2, ncol = 4, byrow = T)
  g <- nrow(gene_pop)
  gene_mat_org <- do.call(rbind, lapply(1:g, function(x){
    pos <- stats::runif(d_each)
    cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(d_each, sd = 0.05),
          pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(d_each, sd = 0.05))
  }))
  d <- nrow(gene_mat_org)

  # form observations
  gram_mat <- cell_mat_org %*% t(gene_mat_org) #natural parameter
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
       h = h, g = g, k = k)
}

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green

#############################

n_seq <- round(seq(10, 50, length.out = 4))
d_seq <- round(seq(20, 100, length.out = 4))
res_list <- vector("list", 4)

for(i in 1:4){
  set.seed(10)
  res <- .data_generator(n_each = n_seq[i], d_each = d_seq[i], multiplier = 0.1)
  dat <- res$dat

  # real analysis
  init <- .initialization(dat, family = "gaussian", max_val = 10)
  res_list[[i]] <- .fit_factorization(dat, init$u_mat, init$v_mat,
                            max_val = 5, family = "gaussian", verbose = T,
                            max_iter = 10, reparameterize = T,
                            return_path = F)
}

png("../figure/simulation/example_trajectories.png", height = 800, width = 3000, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,1))
set.seed(10)
n_pop <- 500
res <- .data_generator(n_each = n_pop, d_each = 50, multiplier = 0.1)
tmp <- res$cell_mat; svd_res <- svd(tmp); tmp <- svd_res$u[,1:res$k] %*% diag(svd_res$d[1:res$k])
plot(tmp[,1], tmp[,2], pch = 16, col = col_vec[rep(1:4, each = n_pop)], asp = T)

for(i in c(1,2,4)){
  tmp <- res_list[[i]]$u_mat; svd_res <- svd(tmp); tmp <- svd_res$u[,1:res$k] %*% diag(svd_res$d[1:res$k])
  plot(tmp[,1], tmp[,2],
       pch = 16, col = col_vec[rep(1:4, each = n_seq[i])], asp = T,
       xlab = "X[,1]", ylab = "X[,2]")
  curves <- slingshot(tmp, cluster_labels = rep(1:4, each = n_seq[i]),
                      starting_cluster = 1)
}
graphics.off()
