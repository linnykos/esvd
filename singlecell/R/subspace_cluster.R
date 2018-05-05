.lasso_connectivity <- function(dat, sparsity = 5){
  n <- nrow(dat)

  sapply(1:n, function(i){
    print(i)
    reg <- glmnet::glmnet(t(dat[-i,,drop = F]), dat[i,], intercept = F)
    sparsity_idx <- which(reg$df <= sparsity)

    vec <- rep(0, n)
    vec[-i] <- as.numeric(glmnet::coef.glmnet(reg, s = reg$lambda[rev(sparsity_idx)[1]])[-1])
    vec
  })
}

.one_dimension_connectivity <- function(dat){
  stopifnot(ncol(dat) == 2)
  n <- nrow(dat)

  l2_vec <- apply(dat, 1, function(x){
    .l2norm(c(x[1], -1/x[2]))
  })

  sapply(1:n, function(i){
    dis_vec <- abs(c(dat[i,1] * dat[-i,1] - dat[i,2] / dat[-i,2]))
    dis_vec <- vec/l2_vec[-i]
    idx <- which.min(dis_vec)

    vec <- rep(0, n)
    vec[-i][idx] <- 1
    vec
  })
}

.coef_to_adj <- function(mat, tol = 1e-6){
  mat <- abs(mat) + t(abs(mat))

  mat[which(abs(mat) >= tol)] <- 1
  mat[which(abs(mat) < tol)] <- 0
  mat
}

.graph_from_adj <- function(adj){
  igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
}

.l2norm <- function(x){
  sqrt(sum(x^2))
}
