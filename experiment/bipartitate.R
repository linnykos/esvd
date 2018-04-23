n = 100
d = 100

prop1 = c(0.5, 0.5)
prop2 = c(0.7, 0.3)
B = matrix(c(0.7, 0.3, 0.3, 0.7), 2, 2)

class_n <- sapply(1:n, function(x){sample(1:2, 1, prob = prop1)})
class_d <- sapply(1:d, function(x){sample(1:2, 1, prob = prop2)})

set.seed(10)
adj <- matrix(0, n, d)
for(i in 1:n){
  for(j in 1:d){
    adj[i,j] = rbinom(1, 1, B[class_n[i], class_d[j]])
  }
}

adj_full <- matrix(0, n+d, n+d)
adj_full[1:n, (n+1):(n+d)] <- adj
adj_full[(n+1):(n+d), 1:n] <- t(adj)

eig <- eigen(adj_full)

plot(sqrt(eig$values[1])*eig$vectors[,1],
     sqrt(abs(eig$values[length(eig$values)]))*eig$vectors[,length(eig$values)],
     asp = T)

col_vec <- c(c("red", "blue")[class_n], c("green", "black")[class_d])

plot(sqrt(eig$values[1])*eig$vectors[,1], sqrt(eig$values[2])*eig$vectors[,2],
     asp = T, col = col_vec, pch = 16)

library(igraph)
g <- igraph::graph_from_adjacency_matrix(adj_full, mode = "undirected")

set.seed(10)
l <- igraph::layout_nicely(g)

igraph::plot.igraph(g, vertex.label = NA, main = "Full graph")
