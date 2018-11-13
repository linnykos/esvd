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

  res <- .identification(cell_mat, gene_mat)
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

set.seed(10)
res <- .data_generator(n_each = 25, d_each = 50, multiplier = 0.1)
stopifnot(sum(abs(t(res$cell_mat)%*%res$cell_mat - t(res$gene_mat)%*%res$gene_mat)) <= 1e-6)
stopifnot(sum(abs((t(res$cell_mat)%*%res$cell_mat)[2:3])) <= 1e-6)
dat <- res$dat

init <- .initialization(dat, family = "gaussian", max_val = 10)
u_mat = init$u_mat
v_mat = init$v_mat
max_val = 5
family = "gaussian"
reparameterize = T
extra_weights = rep(1, nrow(dat))
tol = 1e-3
max_iter = 10
verbose = T
return_path = T
cores = NA

k <- ncol(u_mat)
if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

idx <- which(!is.na(dat))
min_val <- min(dat[which(dat > 0)])
dat[which(dat == 0)] <- min_val/2

current_obj <- Inf
next_obj <- .evaluate_objective(dat, u_mat, v_mat, extra_weights = extra_weights)
obj_vec <- c(next_obj)
if(verbose) print(paste0("Finished initialization : Current objective is ", next_obj))
if(return_path) res_list <- list(list(u_mat = u_mat, v_mat = v_mat)) else res_list <- NA

current_obj <- next_obj

u_mat <- .optimize_mat(dat, u_mat, v_mat, left = T, max_val = max_val, extra_weights = extra_weights,
                       !is.na(cores))
v_mat <- .optimize_mat(dat, v_mat, u_mat, left = F, max_val = max_val, extra_weights = extra_weights,
                       !is.na(cores))

print(range(u_mat %*% t(v_mat)))

if(reparameterize){
  tmp <- .reparameterize(u_mat, v_mat)
  u_mat <- tmp$X; v_mat <- tmp$Y
  print(range(u_mat %*% t(v_mat)))
}

next_obj <- .evaluate_objective(dat, u_mat, v_mat, extra_weights = extra_weights)


