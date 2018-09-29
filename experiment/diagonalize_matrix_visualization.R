rm(list=ls())

construct_A <- function(cov_x, cov_y){
  eigen_x <- eigen(cov_x)
  eigen_y <- eigen(cov_y)

  Ux <- eigen_x$vectors
  Uy <- eigen_y$vectors
  Dx <- diag(eigen_x$values)
  Dy <- diag(eigen_y$values)

  inner_left <- sqrt(Dx) %*% t(Ux) %*% Uy %*% sqrt(Dy)
  svd_res <- svd(inner_left)
  Q <- svd_res$v %*% t(svd_res$u)

  sym_prod <- Ux %*% t(Ux) %*% Uy %*% sqrt(Dy) %*% Q %*% diag(1/sqrt(diag(Dx))) %*% t(Ux)

  eigen_sym <- eigen(sym_prod)
  diag(sqrt(eigen_sym$values)) %*% t(eigen_sym$vectors)
}

dmvnorm <- function(x, mean, sigma){
  as.numeric(exp(-.5*t(x-mean)%*%solve(sigma)%*%(x-mean)) /
               sqrt((2*pi)^length(mean)*det(sigma)))
}

#create a density over a grid
create_prob_grid <- function(xlim, ylim, density_func, grid_size = 100){
  xseq <- seq(xlim[1], xlim[2], length.out = grid_size)
  yseq <- seq(ylim[1], ylim[2], length.out = grid_size)

  mat <- sapply(xseq, function(x){ sapply(yseq,
                                          function(y){density_func(x,y)})})

  rownames(mat) <- yseq; colnames(mat) <- xseq
  mat <- mat[nrow(mat):1,]
  mat/sum(mat)

  list(mat = mat, func = density_func)
}

rotate <- function(mat){t(mat)[,nrow(mat):1]}
l2norm <- function(x){sqrt(sum(x^2))}

#plot the density
plot_grid <- function(obj, type = c("contour"), alpha = 0.5, nlevels = 10,
                      scaling = 1, ...){
  mat <- obj$mat
  xseq <- as.numeric(colnames(mat)); yseq <- rev(as.numeric(rownames(mat)))

  image(xseq, yseq, t(mat[nrow(mat):1,]), col = heat.colors(100, alpha = alpha),
        asp = 1, xlab = "X", ylab = "Y", ...)
  contour(xseq, yseq, t(mat[nrow(mat):1,]), nlevels = nlevels, add = T, drawlabels = F, ...)

  lines(rep(0,2), range(yseq), lwd = 2, lty = 2)
  lines(range(xseq), rep(0,2), lwd = 2, lty = 2)

  sm <- compute_second_moment(obj)
  eig <- eigen(sm)

  vec_list <- lapply(1:2, function(i){
    scaling*eig$values[i]*eig$vectors[,i]
  })

  for(i in 1:2){
    arrows(0,0, vec_list[[i]][1], vec_list[[i]][2], col = "blue", length = 0.1, lwd = 3, angle = 10, ...)
    arrows(0,0, -vec_list[[i]][1], -vec_list[[i]][2], col = "blue", length = 0.1, lwd = 3, angle = 10, ...)
  }

  invisible()
}

compute_mean <- function(obj){
  mat <- obj$mat
  col_vec <- colSums(mat)
  val1 <- as.numeric(colnames(mat)) %*% col_vec / sum(col_vec)

  row_vec <- rowSums(mat)
  val2 <- as.numeric(rownames(mat)) %*% row_vec / sum(row_vec)

  c(val1, val2)
}

compute_covariance <- function(obj){
  mat <- obj$mat
  mean_vec <- compute_mean(obj)

  x_centered <- as.numeric(colnames(mat)) - mean_vec[1]
  y_centered <- as.numeric(rownames(mat)) - mean_vec[2]
  prod_mat <- x_centered %*% t(y_centered)
  diag_val <- sum(mat * prod_mat)/sum(mat)

  x_weight <- colSums(mat)
  x_val <- x_weight %*% x_centered^2 / sum(x_weight)

  y_weight <- rowSums(mat)
  y_val <- y_weight %*% y_centered^2 / sum(y_weight)

  cov_mat <- matrix(diag_val, 2, 2)
  cov_mat[1,1] <- x_val; cov_mat[2,2] <- y_val

  cov_mat
}

compute_second_moment <- function(obj){
  mat <- obj$mat
  mean_vec <- compute_mean(obj)

  x_seq <- as.numeric(colnames(mat))
  y_seq <- as.numeric(rownames(mat))
  prod_mat <- x_seq %*% t(y_seq)
  diag_val <- sum(mat * prod_mat)/sum(mat)

  x_weight <- colSums(mat)
  x_val <- x_weight %*% x_seq^2 / sum(x_weight)

  y_weight <- rowSums(mat)
  y_val <- y_weight %*% y_seq^2 / sum(y_weight)

  cov_mat <- matrix(diag_val, 2, 2)
  cov_mat[1,1] <- x_val; cov_mat[2,2] <- y_val

  cov_mat
}

apply_transformation <- function(obj, A){
  mat <- obj$mat
  xseq <- as.numeric(colnames(mat)); yseq <- as.numeric(rownames(mat))
  len_x <- length(xseq); len_y <- length(yseq)

  elementwise <- cbind(rep(xseq, each = len_y), rep(yseq, times = len_x))
  preimage <- elementwise %*% solve(t(A))

  prob_mass <- apply(preimage, 1, function(x){
    obj$func(x[1], x[2])
  })

  new_mat <- matrix(prob_mass, ncol = ncol(mat), nrow = nrow(mat))
  colnames(new_mat) <- colnames(mat)
  rownames(new_mat) <- rownames(mat)

  new_func <- function(x, y){
    val <- solve(A) %*% c(x,y)
    obj$func(val[1], val[2])
  }

  list(mat = new_mat, func = new_func)
}

##################

set.seed(10)

density_func1 <- function(x,y){
  val1 <- dmvnorm(c(x,y), mean = c(-1.5,1.5), sigma = matrix(c(1,-.5,-.5,1),2,2))

  val1
}

density_func2 <- function(x,y){
  val1 <- dmvnorm(c(x,y), mean = c(1,0.5), sigma = matrix(c(.5,0,0,3),2,2))

  val1
}

xlim <- c(-5, 5); ylim <- c(-5, 5)
obj1 <- create_prob_grid(xlim, ylim, density_func1, grid_size = 100)
obj2 <- create_prob_grid(xlim, ylim, density_func2, grid_size = 100)

par(mfrow = c(1,2))
plot_grid(obj1, main = "Distribution 1, before")
plot_grid(obj2, main = "Distribution 2, before")

sm1 <- compute_second_moment(obj1)
sm2 <- compute_second_moment(obj2)

A <- construct_A(sm1, sm2)
new_obj1 <- apply_transformation(obj1, A)
new_obj2 <- apply_transformation(obj2, solve(t(A)))

par(mfrow = c(1,2))
plot_grid(new_obj1, main = "Distribution 1, after")
plot_grid(new_obj2, main = "Distribution 2, after")

new_sm1 <- compute_second_moment(new_obj1)
new_sm2 <- compute_second_moment(new_obj2)
new_sm1; new_sm2
