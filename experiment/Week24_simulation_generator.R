.data_generator <- function(distr_func = function(x){stats::rnorm(1, 1/x, 2/x)},
                            min_val = 1e-4,
                            total = 150, col_drop = T,
                            n_each = 50){

  #construct the cell information
  cell_pop <- matrix(c(4,10, 25,100,
                       40,10, 60,80,
                       60,80, 25,100,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  cell_mat_org <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  n <- nrow(cell_mat_org)
  k <- ncol(cell_mat_org)

  # construct the gene information
  gene_pop <- matrix(c(20, 90, 25, 100,
                          90,20, 100,25)/10, nrow = 2, ncol = 4, byrow = T)
  g <- nrow(gene_pop)
  d_each <- 120
  gene_mat_org <- do.call(rbind, lapply(1:g, function(x){
    pos <- stats::runif(d_each)
    cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(d_each, sd = 0.1),
          pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(d_each, sd = 0.1))
  }))
  d <- nrow(gene_mat_org)

  # form observations
  gram_mat <- cell_mat_org %*% t(gene_mat_org) #natural parameter
  inv_mat <- 1/gram_mat #the actual mean
  inv_mat[inv_mat >= 1] <- 1
  inv_mat[inv_mat <= 1e-6] <- 1e-6
  gram_mat <- 1/inv_mat
  svd_res <- svd(gram_mat)
  cell_mat <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
  gene_mat <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
  gram_mat <- cell_mat %*% t(gene_mat)

  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- distr_func(max(gram_mat[i,j], 1e-4))
    }
  }

  # now start fidgeting with dropout
  obs_mat2 <- obs_mat
  threshold <- quantile(as.numeric(obs_mat[obs_mat > 0]), probs = 0.05)
  obs_mat2[obs_mat2 < threshold] <- 0 # the "true zeros"
  length(which(obs_mat2 == 0))/prod(dim(obs_mat2))

  # now do something more dramatic with dropout
  obs_mat3 <- obs_mat2
  .dropped_indices <- function(x, total){
    vec <- 1:length(x)
    samp <- sample(vec, size = total, replace = T, prob = x)
    setdiff(vec, unique(samp))
  }

  if(col_drop){
    total_vec <- rep(total, ncol(obs_mat3))
    for(i in 1:ncol(obs_mat3)){
      idx <- .dropped_indices(obs_mat3[,i], total = total_vec[i])
      obs_mat3[idx,i] <- 0
    }
  } else {
    total_vec <- rep(total, nrow(obs_mat3))
    for(i in 1:nrow(obs_mat3)){
      idx <- .dropped_indices(obs_mat3[i,], total = total_vec[i])
      obs_mat3[i,idx] <- 0
    }
  }
  length(which(obs_mat3 == 0))/prod(dim(obs_mat3))

  list(dat = obs_mat3,
       cell_pop = cell_pop, cell_mat = cell_mat, gene_mat = gene_mat,
       cell_mat_org = cell_mat_org, gene_mat_org = gene_mat_org,
       gram_mat = gram_mat, obs_mat = obs_mat, obs_mat2 = obs_mat2,
       n_each = n_each, d_each = d_each, h = h, g = g, k = k)
}
