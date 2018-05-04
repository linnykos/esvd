.lasso_connectivity <- function(dat, sparsity = 5){
  n <- nrow(dat)

  coef_mat <- sapply(1:n, function(i){
    print(i)
    reg <- glmnet::glmnet(t(dat[-i,,drop = F]), dat[i,], intercept = F)
    sparsity_idx <- which(reg$df <= sparsity)

    vec <- rep(0, n)
    vec[-i] <- as.numeric(glmnet::coef.glmnet(reg, s = reg$lambda[rev(sparsity_idx)[1]])[-1])
    vec
  })

  abs(coef_mat) + t(abs(coef_mat))
}

.coef_to_adj <- function(mat, tol = 1e-6){
  mat[which(abs(mat) >= tol)] <- 1
  mat[which(abs(mat) < tol)] <- 0
  mat
}

.graph_from_adj <- function(adj){
  igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
}
