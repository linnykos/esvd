rm(list=ls())

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

generate_latent <- function(res, proportion, n){
  lapply(1:n, function(x){
    idx <- sample(1:length(proportion), 1, prob = proportion)
    vec <- list(positive = res$positive[,idx], negative = res$negative[,idx],
                class = idx)
  })
}

inner_prod <- function(vec1, vec2){
  as.numeric(sqrt((vec1$positive%*%vec2$positive) + (vec1$negative%*%vec2$negative)))
}

generate_graph <- function(dat){
  n <- length(dat)
  adj <- matrix(0, n, n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      adj[i,j] <- rbinom(1, 1, inner_prod(dat[[i]], dat[[j]]))
      adj[j,i] <- adj[i,j]
    }
  }

  adj
}

##########
proportion <- c(0.3, 0.3, 0.4)
# B <- matrix(0, 3, 3)
# for(i in 1:3){
#   for(j in 1:3){
#     B[i,j] <- .25 + .25*(i != j)
#   }
# }
B <- matrix(c(1/4, 1/2, 1/4, 1/2, 1/4, 1/4, 1/4, 1/4, 1/6), 3, 3)
#B <- matrix(c(3/4, 1/4, 1/4, 1/4, 3/4, 1/4, 1/4, 1/4, 3/4), 3, 3)

set.seed(10)
res <- decompose(B)
dat <- generate_latent(res, proportion, 200)
adj <- generate_graph(dat)

library(igraph)
g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")

set.seed(10)
l <- igraph::layout_nicely(g)
node_col <- sapply(dat, function(x){c("red","blue","green")[x$class]})

igraph::plot.igraph(g, layout = l, vertex.label = NA,
                    vertex.color = node_col, main = "Full graph")

xlim <- range(sapply(dat, function(x){x$positive}))
ylim <- range(sapply(dat, function(x){x$negative}))

plot(NA, xlim = xlim, ylim = ylim)
for(i in 1:length(dat)){
  points(dat[[i]]$positive, dat[[i]]$negative, pch = 16)
}

