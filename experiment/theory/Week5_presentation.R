rm(list=ls())
lambda_func <- function(x,y){
  val1 <- 0.5*mvtnorm::dmvnorm(c(x,y), mean = c(.1, .1), sigma = matrix(c(.5,.25,.25,.25),2,2)/5)
  val2 <- 0.5*mvtnorm::dmvnorm(c(x,y), mean = c(0.75,0.9), sigma = matrix(c(0.25,0,0,.25),2,2)/10)

  val3 <- 0.5*mvtnorm::dmvnorm(c(y,x), mean = c(.1, .1), sigma = matrix(c(.5,.25,.25,.25),2,2)/5)
  val4 <- 0.5*mvtnorm::dmvnorm(c(y,x), mean = c(0.75,0.9), sigma = matrix(c(.25,0,0,.25),2,2)/10)

  val1 + val2 + val3 + val4
}

create_prob_grid <- function(density_func, xlim = c(0,1), ylim = c(0,1), grid_size = 100){
  xseq <- seq(xlim[1], xlim[2], length.out = grid_size)
  yseq <- seq(ylim[1], ylim[2], length.out = grid_size)

  mat <- sapply(xseq, function(x){ sapply(yseq,
                                          function(y){density_func(x,y)})})

  rownames(mat) <- yseq; colnames(mat) <- xseq
  mat <- mat[nrow(mat):1,]
  mat/sum(mat)

  list(mat = mat, func = density_func)
}

plot_grid <- function(mat, type = c("contour"), alpha = 0.5, nlevels = 10,
                      scaling = 1, ...){
  xseq <- as.numeric(colnames(mat)); yseq <- rev(as.numeric(rownames(mat)))

  graphics::image(xseq, yseq, t(mat[nrow(mat):1,]), col = rev(gray.colors(100, alpha = alpha)),
        asp = 1, ...)
  contour(xseq, yseq, t(mat[nrow(mat):1,]), nlevels = nlevels, add = T, drawlabels = F)

  invisible()
}


rotate <- function(mat){t(mat)[,nrow(mat):1]}

set.seed(10)
grid <- create_prob_grid(lambda_func)
range(grid$mat) #good. just remember to invert it

png("../figure/theory/5_lambda_func.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(2,2,0.5,0.5))
plot_grid(grid$mat, alpha = 1, xlab = "", ylab = "", axes = F)
graphics.off()

#################

# start with the simple graph
generate_graphon <- function(grid, threshold){
  mat <- matrix(0, nrow = nrow(grid$mat), ncol = ncol(grid$mat))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      mat[i,j] <- stats::pexp(threshold, rate = grid$mat[i,j])
    }
  }

  rownames(mat) <- rownames(grid$mat); colnames(mat) <- colnames(grid$mat)

  mat
}

graphon1 <- generate_graphon(grid, 1)
png("../figure/theory/5_graphon1.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(2,2,0.5,0.5))
plot_grid(graphon1, alpha = 1, xlab = "", ylab = "", axes = F, breaks = seq(0,1,length.out = 101))
graphics.off()

graphon0.25 <- generate_graphon(grid, 0.25)
png("../figure/theory/5_graphon0.25.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(2,2,0.5,0.5))
plot_grid(graphon0.25, alpha = 1, xlab = "", ylab = "", axes = F, breaks = seq(0,1,length.out = 101))
graphics.off()

graphon0.5 <- generate_graphon(grid, 0.5)
png("../figure/theory/5_graphon0.5.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(2,2,0.5,0.5))
plot_grid(graphon0.5, alpha = 1, xlab = "", ylab = "", axes = F, breaks = seq(0,1,length.out = 101))
graphics.off()

graphon5 <- generate_graphon(grid, 5)
png("../figure/theory/5_graphon5.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(2,2,0.5,0.5))
plot_grid(graphon5, alpha = 1, xlab = "", ylab = "", axes = F, breaks = seq(0,1,length.out = 101))
graphics.off()

lookup_graphon <- function(graphon, x, y){
  xseq <- colnames(graphon); yseq <- rownames(graphon)

  idx_x <- which.min(abs(x - as.numeric(xseq)))
  idx_y <- which.min(abs(y - as.numeric(yseq)))

  graphon[idx_x, idx_y]
}

# generate 6 points
set.seed(10)
n <- 6
unif_i <- stats::runif(n)
combn_mat <- utils::combn(n, 2)
unif_ij <- stats::runif(ncol(combn_mat))
obs <- matrix(0, n, n)
for(i in 1:(n-1)){
  for(j in (i+1):n){
    threshold <- lookup_graphon(graphon0.5, unif_i[i], unif_i[j])
    val <- unif_ij[intersect(which(combn_mat[1,]==i),which(combn_mat[2,]==j))]
    obs[i,j] <- ifelse(threshold >= val, 1, 0)
    obs[j,i] <- obs[i,j]
  }
}

# plot the three graphs
set.seed(10)
g <- igraph::graph_from_adjacency_matrix(obs, mode = "undirected")
l <- igraph::layout_nicely(g)
l[,1] <- (l[,1]-min(l[,1]))/diff(range(l[,1]))
l[,2] <- (l[,2]-min(l[,2]))/diff(range(l[,2]))

# DOESN'T WORKS
for(i in c(3,4,6)){
  png(paste0("../figure/theory/5_graph_",i,".png"), height = 500, width = 500, res = 300, units = "px")
  par(mar = rep(0.5,4))
  g <- igraph::graph_from_adjacency_matrix(obs[1:i, 1:i], mode = "undirected")
  igraph::plot.igraph(g, layout = l[1:i,], vertex.color = "gray",
                      xlim = c(-1,1), ylim = c(-1,1))
  graphics.off()
}

real_obs <- matrix(0, n, n)
for(i in 1:(n-1)){
  for(j in (i+1):n){
    lambda <- lookup_graphon(grid$mat, unif_i[i], unif_i[j])
    val <- unif_ij[intersect(which(combn_mat[1,]==i),which(combn_mat[2,]==j))]
    real_obs[i,j] <- stats::qexp(val, rate = lambda)
    real_obs[j,i] <- real_obs[i,j]
  }
}

#####################

lambda <- lookup_graphon(grid$mat, 0.1, 0.4)
x_seq <- seq(0, 2, length.out = 100)
y_seq <- sapply(x_seq, function(x){stats::pexp(x, rate = lambda)})

png("../figure/theory/5_exponential.png", height = 1000, width = 1200, res = 300, units = "px")
plot(x_seq, y_seq, type = "l", axes = F, xlab = "Threshold", lwd = 2,
     ylab = "Graphon value")
axis(1)
axis(2)
for(i in c(0.25, 0.5, 1, 1.5)){
  idx <- which.min(abs(x_seq - i))
  lines(rep(x_seq[idx], 2), c(0, y_seq[idx]), col = rgb(205,40,54,maxColorValue=255),
        lty = 2, lwd = 2)
  points(x_seq[idx], y_seq[idx], col = rgb(205,40,54,maxColorValue=255), pch = 16, cex = 1)
}
graphics.off()
