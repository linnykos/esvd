rm(list=ls())

# the function to generate data
generator <- function(n_each, d_each){
  sigma <- 0.01

  cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                       40,10, 60,80, 60,80, 100, 25)/10,
                     nrow = 4, ncol = 4, byrow = T)
  gene_pop <- matrix(c(20,90, 25,100,
                       90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)

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

  alpha <- 2
  obs_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- rnorm(1, 1/nat_mat[i,j], sd = 1/(alpha*nat_mat[i,j]))
    }
  }

  list(obs_mat = obs_mat, nat_mat = nat_mat,
       cell_mat = cell_mat, gene_mat = gene_mat)
}


####################

# function for objective, gradient and heissan

heissan <- function(nat_mat, x_mat, y_mat, alpha = 2){
  stopifnot(ncol(x_mat) == ncol(y_mat))
  n <- nrow(x_mat)
  p <- nrow(y_mat)
  k <- ncol(x_mat)

  a_mat <- (1/alpha^2 + 1)*1/(nat_mat^2)

  heis <- matrix(0, n*k, n*k)
  for(i in 1:n){
    idx <- ((i-1)*k+1):(i*k)
    for(j in 1:p){
      # tmp <- (1/(n*p) * as.numeric((1/(x_mat[i,] %*% y_mat[j,])) +
      #                      alpha^2*a_mat[i,j]) * (y_mat[j,] %*% t(y_mat[j,])))
      tmp <- (1/(p) * as.numeric((1/(x_mat[i,] %*% y_mat[j,])) +
                                     alpha^2*a_mat[i,j]) * (y_mat[j,] %*% t(y_mat[j,])))
      heis[idx, idx] <- tmp + heis[idx, idx]
    }
  }

  heis
}

###############

# res <- generator(50, 50)
# zz <- heissan(res$nat_mat, res$cell_mat, res$gene_mat)
# yy <- eigen(zz)$values

###################

zz <- sapply(round(seq(10,500, length.out = 11)), function(x){
  print(x)
  res <- generator(x, x)
  heis <- heissan(res$nat_mat, res$cell_mat, res$gene_mat)
  range(eigen(heis)$values)
})

zz
