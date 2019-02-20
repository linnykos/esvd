rm(list=ls())
## from http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

# visualize a 2-dimensional normal
ellipse_points <- function(mean_vec, cov_mat){
  eig <- eigen(cov_mat)
  alpha <- atan(eig$vectors[2,1]/eig$vectors[1,1])
  if(alpha < 0) alpha <- alpha + 2*pi

  a <- sqrt(eig$values[1])
  b <- sqrt(eig$values[2])

  theta_grid <- seq(0, 2*pi, length.out = 100)
  ellipse_x <- a*cos(theta_grid)
  ellipse_y <- b*sin(theta_grid)

  R <- matrix(c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)), 2, 2)
  val <- cbind(ellipse_x, ellipse_y) %*% R

  val <- t(apply(val, 1, function(x){x + mean_vec}))
}

dist_func <- function(mean_vec1, cov_mat1, mean_vec2, cov_mat2){
  as.numeric(t(mean_vec1 - mean_vec2) %*% solve(cov_mat1 + cov_mat2) %*% (mean_vec1 - mean_vec2))
}

visualization <- function(mean_vec1, cov_mat1, mean_vec2, cov_mat2){

  set.seed(10)

  e1 <- ellipse_points(mean_vec1, cov_mat1)
  x1 <- MASS::mvrnorm(100, mean_vec1, cov_mat1)
  e2 <- ellipse_points(mean_vec2, cov_mat2)
  x2 <- MASS::mvrnorm(100, mean_vec2, cov_mat2)

  plot(x1[,1], x1[,2], pch = 16, asp = T, xlim = range(c(x1[,1], x2[,1])),
       ylim = range(c(x1[,2], x2[,2])), col = rgb(0,0,0,0.2),
       main = paste0(dist_func(mean_vec1, cov_mat1, mean_vec2, cov_mat2)))
  points(x2[,1], x2[,2], pch = 16, col = rgb(1,0,0,0.2))

  points(e1[,1], e1[,2], cex = 0.5, col = "black")
  points(e2[,1], e2[,2], cex = 0.5, col = "red")

}

mean_vec1 <- c(2,-1)
cov_mat1 <- matrix(c(2,-1.9,-1.9,2),2,2)
mean_vec2 <- c(1,1)
cov_mat2 <- matrix(c(2,-1.9,-1.9,2),2,2)
visualization(mean_vec1, cov_mat1, mean_vec2, cov_mat2)
