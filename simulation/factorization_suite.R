rm(list=ls())
library(simulation)
library(singlecell)

paramMat <- cbind(round(exp(seq(log(10), log(200), length.out = 10))),
                  round(exp(seq(log(20), log(400), length.out = 10))))
colnames(paramMat) <- c("n", "d")
trials <- 50

################

# setup
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100,25)/50,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/50, nrow = 2, ncol = 4, byrow = T)
distr_func = function(x){stats::rnorm(1, 4/x, sd = 2/x)}
n_each = 50
d_each = 100
sigma = 0.05
total = 150

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

obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
for(i in 1:n){
  for(j in 1:d){
    obs_mat[i,j] <- distr_func(max(gram_mat[i,j], 1e-4))
  }
}

obs_mat[obs_mat < 0] <- 0
obs_mat2 <- obs_mat

# now do something more dramatic with dropout
obs_mat3 <- obs_mat2
.dropped_indices <- function(x, total){
  vec <- 1:length(x)
  samp <- sample(vec, size = total, replace = T, prob = x)
  setdiff(vec, unique(samp))
}

total_vec <- rep(total, nrow(obs_mat3))
for(i in 1:nrow(obs_mat3)){
  idx <- .dropped_indices(obs_mat3[i,], total = total_vec[i])
  obs_mat3[i,idx] <- 0
}
length(which(obs_mat3 == 0))/prod(dim(obs_mat3))
quantile(obs_mat3)

######xx

obs_mat4 <- pmin(obs_mat3, 10)
obs_mat4 <- round(exp(obs_mat4))-1
quantile(obs_mat4)
quantile(obs_mat4[obs_mat4>0])

impute_res <- SAVER::saver(t(obs_mat4))
zz <- t(impute_res$estimate)
zz <- log(zz+1)
plot(as.numeric(zz), as.numeric(obs_mat2), asp = T)

