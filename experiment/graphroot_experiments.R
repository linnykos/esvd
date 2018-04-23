B = matrix(0, 3, 3)
for(i in 1:3){
  for(j in 1:3){
    B[i,j] = .25 + .25*(i != j)
  }
}

xx = c(0.65, 0.41, 0)
yy = c(0.65, -0.2, -.35)
zz = c(0.65, -.2, 0.35)

innerprod <- function(vec1, vec2){
  sqrt((vec1[1]*vec2[1]) + (vec1[2:3]%*%vec2[2:3]))
}

#how did that happen??

B <- matrix(c(1/4, 1/2, 1/4, 1/2, 1/4, 1/4, 1/4, 1/4, 1/6), 3, 3)
marginal <- c(0.3, 0.3, 0.4)

generate_adj <- function(B, marginal, n = 100){
  idx <- rmultinom(n, 1, marginal)
  membership <- t(idx)%*%B%*%idx
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      membership[i,j] <- stats::rbinom(1, 1, membership[i,j])
      membership[j,i] <- membership[i,j]
    }
  }

  diag(membership) <- 0
  list(adj = membership, label = apply(idx, 2, function(x){which(x != 0)}))
}

res <- generate_adj(B, marginal, n = 1000)
eig <- eigen(res$adj)

plot(sqrt(eig$values[1])*eig$vectors[,1],
     sqrt(abs(eig$values[length(eig$values)]))*eig$vectors[,length(eig$values)],
     asp = T)

decompose <- function(B, tol = 1e-6){
  eig <- eigen(B)
  eig_val <- sqrt(abs(eig$values))
  pos_idx <- which(eig$values >= tol)
  neg_idx <- which(eig$values <= -tol)

  lis <- list(positive = matrix(0, length(pos_idx), nrow(B)),
              negative = matrix(0, length(neg_idx), nrow(B)))

  for(i in 1:nrow(B)) {
    lis$positive[,i] <- eig_val[pos_idx] * eig$vectors[i,pos_idx]
    lis$negative[,i] <- eig_val[neg_idx] * eig$vectors[i,neg_idx]
  }

  lis
}

library(igraph)
g <- igraph::graph_from_adjacency_matrix(res$adj, mode = "undirected")

set.seed(10)
l <- igraph::layout_nicely(g)
node_col <- c("red","blue","green")[res$label]

igraph::plot.igraph(g, vertex.label = NA,
                    vertex.color = node_col, main = "Full graph")
