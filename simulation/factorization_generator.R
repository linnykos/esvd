generate_natrual_mat <- function(cell_pop, gene_pop, n_each, d_each, sigma){
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

  nat_mat <- cell_mat %*% t(gene_mat)

  list(nat_mat = nat_mat, cell_mat = cell_mat, gene_mat = gene_mat)
}

generate_dropout <- function(obs_mat, total){
  obs_mat3 <- obs_mat
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

  obs_mat3
}

#####################

generator_exponential <- function(nat_mat){
  obs_mat <- matrix(0, ncol = ncol(nat_mat), nrow = nrow(nat_mat))
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- rexp(1, nat_mat[i,j])
    }
  }

  obs_mat[obs_mat < 0] <- 0

  obs_mat
}


