.data_generator1 <- function(cell_addition, gene_addition,
                            distr_func = function(x){stats::rnorm(1, x, sd = 1/2)},
                            rotate_cells = F){

  #construct the cell information
  cell_pop <- matrix(c(0, -1/4, -1/2, -5/4,
                       0, -1/4, 1/2, -5/4,
                       0, -1/4, 0, 0,
                       0, 0, 0, 1), nrow = 4, ncol = 4, byrow = T) + cell_addition
  h <- nrow(cell_pop)
  n_each <- 50
  cell_mat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 1/10),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 1/10))
  }))
  n <- nrow(cell_mat)
  k <- ncol(cell_mat)

  #now diagonalize it
  if(rotate_cells){
    cell_mat <- cell_mat %*% svd(cell_mat)$v
    cell_mat[,1] <- cell_mat[,1] - min(cell_mat[,1])
    cell_mat[,2] <- cell_mat[,2] - min(cell_mat[,2])
  }
  # check: t(cell_mat) %*% cell_mat / n

  # construct the gene information
  g <- 2
  d_each <- 120
  gene_pop <- diag(2)

  gene_mat <- do.call(rbind, lapply(1:g, function(x){
    MASS::mvrnorm(n = d_each, mu = gene_pop[x,], Sigma = 0.25*diag(2))
  }))
  gene_mat <- gene_mat + gene_addition
  d <- nrow(gene_mat)

  # gene_mat <- gene_mat %*% svd(gene_mat)$v

  # form observations
  gram_mat <- cell_mat %*% t(gene_mat)
  obs_mat <- matrix(0, ncol = ncol(gram_mat), nrow = nrow(gram_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- distr_func(gram_mat[i,j])
    }
  }

  obs_mat2 <- obs_mat

  # now start fidgeting with dropout
  obs_mat2[obs_mat2 < 0] <- 0

  # now do something more dramatic with dropout
  .dropped_indices <- function(x, total){
    vec <- 1:length(x)
    samp <- sample(vec, size = total, replace = T, prob = x)
    setdiff(vec, unique(samp))
  }

  if(rotate_cells){
    total_vec <- round(sample(seq(50, 300, length.out = nrow(obs_mat2))))
  } else {
    total_vec <- rep(200, nrow(obs_mat2))
  }

  for(i in 1:nrow(obs_mat2)){
    idx <- .dropped_indices(obs_mat2[i,], total = total_vec[i])
    obs_mat2[i,idx] <- 0
  }
  length(which(obs_mat2 == 0))/prod(dim(obs_mat2))

  list(dat = obs_mat2,
       cell_pop = cell_pop, cell_mat = cell_mat, gene_mat = gene_mat,
       gram_mat = gram_mat, obs_mat = obs_mat,
       n_each = n_each, d_each = d_each, h = h, g = g, k = k)
}
