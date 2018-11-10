rm(list=ls())
library(singlecell)

.data_generator <- function(distr_func = function(x){stats::rnorm(1, 1/x, sd = sqrt(1/(2*x)))},
                            n_each = 50, d_each = 120, sigma = 0.05,
                            multiplier = 1){

  #construct the cell information
  cell_pop <- multiplier*matrix(c(4,10, 25,100,
                       40,10, 60,80,
                       60,80, 25,100,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
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

set.seed(10)
res <- .data_generator(n_each = 20, d_each = 40, multiplier = 0.1)
dat <- res$dat
quantile(dat)
# .plot_singlecell(res$dat)
plot(res$cell_mat[,1], res$cell_mat[,2], pch = 16, asp = T,
     col = col_vec[rep(1:4, each = res$n_each)])

# naive analysis
svd_res <- svd(dat)
u_mat <- svd_res$u[,1:res$k] %*% diag(sqrt(svd_res$d[1:res$k]))
plot(u_mat[,1], u_mat[,2], pch = 16, asp = T,
     col = col_vec[rep(1:4, each = res$n_each)])

# real analysis
init <- .initialization(dat, family = "gaussian", max_val = 10)
plot(init$u_mat[,1], init$u_mat[,2], pch = 16, asp = T,
     col = col_vec[rep(1:4, each = res$n_each)])
fit_a <- .fit_factorization(dat, init$u_mat, init$v_mat,
                          max_val = 5, family = "gaussian", verbose = T,
                          max_iter = 50)
plot(fit_a$u_mat[,1], fit_a$u_mat[,2], pch = 16, asp = T,
     col = col_vec[rep(1:4, each = res$n_each)])
fit_b <- .fit_factorization(dat, init$u_mat, init$v_mat,
                          max_val = 5, family = "gaussian", verbose = T,
                          max_iter = 50, reparameterize = T)
plot(fit_b$u_mat[,1], fit_b$u_mat[,2], pch = 16, asp = T,
     col = col_vec[rep(1:4, each = res$n_each)])

# ideal analysis
fit2 <- .fit_factorization(dat, res$cell_mat, res$gene_mat,
                          max_val = 5, family = "gaussian", verbose = T,
                          max_iter = 50)
plot(fit2$u_mat[,1], fit2$u_mat[,2], pch = 16, asp = T,
     col = col_vec[rep(1:4, each = res$n_each)])
fit2_b <- .fit_factorization(dat, res$cell_mat, res$gene_mat,
                           max_val = 5, family = "gaussian", verbose = T,
                           max_iter = 50, reparameterize = T)
plot(fit2_b$u_mat[,1], fit2_b$u_mat[,2], pch = 16, asp = T,
     col = col_vec[rep(1:4, each = res$n_each)])
