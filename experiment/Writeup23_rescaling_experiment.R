rm(list=ls())
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)
n_each <- 50
d_each <- 200
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

#########

# start the rescaling via linear program
obj <- rep(1, d)
mat <- matrix(0, n-1, d)
for(i in 1:nrow(mat)){
  mat[i,] <- obs_mat[i,] - obs_mat[i+1,]
}
constr_ub <- rep(100, n-1)
constr_lb <- rep(-100, n-1)
var_lb <- rep(0, d)
var_ub <- rep(Inf, d)
zz <- clplite::clp_solve(obj, mat, constr_ub = constr_ub, constr_lb = constr_lb,
                         var_lb = var_lb, var_ub = var_ub, max = TRUE)
zz$solution

