rm(list=ls())
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)
n_each <- 50
d_each <- 100
sigma <- 0.05
total <- 500

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

pred_mat <- cell_mat %*% t(gene_mat)

obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
for(i in 1:n){
  for(j in 1:d){
    obs_mat[i,j] <- rexp(1, pred_mat[i,j])
  }
}

obs_mat[obs_mat < 0] <- 0
obs_mat2 <- round(100*obs_mat)
# obs_mat2 <- round(exp(obs_mat*10)-1)
# length(which(obs_mat2 > 5000))/prod(dim(obs_mat2))
# obs_mat2[obs_mat2 > 200] <- 200
# quantile(obs_mat2, probs = seq(0,1,length.out=11))
# length(which(obs_mat2 == 0))/prod(dim(obs_mat2))

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
# quantile(obs_mat3, probs = seq(0,1,length.out=11))
# length(which(obs_mat3 == 0))/prod(dim(obs_mat3))

plot(cell_mat[,1], cell_mat[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])


##############################

# dat <- obs_mat3
# reweight_factor <- rowSums(dat)
# dat_rescale <- t(sapply(1:nrow(dat), function(i){dat[i,]/reweight_factor[i]}))
# dat_rescale <- dat_rescale*10/mean(dat_rescale)

set.seed(10)
dat_rescale <- obs_mat3
library(NMF)
init <- initialization(dat_rescale, family = "gaussian", k = 2, max_val = -2000)
plot(init$u_mat[,1], init$u_mat[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = n_each)])
res_our <- fit_factorization(dat_rescale, u_mat = init$u_mat, v_mat = init$v_mat,
                                         family = "gaussian",  reparameterize = T,
                                         max_iter = 100, max_val = -2000,
                                         scalar = 2,
                                         return_path = F, cores = 1,
                                         verbose = T)
plot(res_our$u_mat[,1], res_our$u_mat[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

#######

set.seed(10)
zz <- Rtsne::Rtsne(dat_rescale)
plot(zz$Y[,1], zz$Y[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

#####
set.seed(10)
tmp <- vamf:::vamf(t(dat_rescale), 2, log2trans=F) # error in the max-elbos

set.seed(10)
dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_rescale)))
tmp <- zinbwave::zinbwave(dat_se, K = 2, maxiter.optimize = 100)
res_zinb <- tmp@reducedDims$zinbwave

plot(res_zinb[,1], res_zinb[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

########
set.seed(10)
res <- pCMF::pCMF(dat_rescale, K = 2)
plot(res$factor$U[,1], res$factor$U[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = n_each)])


