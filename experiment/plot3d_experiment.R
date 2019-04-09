library(plot3D)
par(mar = c(0, 0, 0, 0))
# Shape 1
M <- mesh(seq(0, 6*pi, length.out = 80),
          seq(pi/3, pi, length.out = 80))
u <- M$x ; v <- M$y
x <- u/2 * sin(v) * cos(u)
y <- u/2 * sin(v) * sin(u)
z <- u/2 * cos(v)
surf3D(x, y, z, colvar = z, colkey = FALSE)

surf3D(x[,1:2], y[,1:2], z[,1:2], colvar = z[,1:2], colkey = FALSE)
surf3D(x[,1:2], y[,1:2], z[,1:2], colvar = z[,1:2], colkey = FALSE, theta = 0, phi = 180)
idx <- sample(1:nrow(x))
surf3D(x[idx,1:2], y[idx,1:2], z[idx,1:2], colvar = z[idx,1:2], colkey = FALSE)

#############

# let's try to make a straight tube?
seq_vec <- seq(0, 2*pi, length.out = 20)
x <- cos(seq_vec); y <- sin(seq_vec)
x_mat <- do.call(cbind, lapply(1:10, function(i){x}))
y_mat <- do.call(cbind, lapply(1:10, function(i){y}))
z_mat <- do.call(rbind, lapply(1:20, function(i){1:10}))

surf3D(x_mat, y_mat, z_mat, colvar = z_mat, bty = "f")

#############
#
# # now let's try a curved tube
seq_vec <- seq(0, 4*pi, length.out = 100)
x <- seq_vec; y <- 5*cos(seq_vec)
dat <- cbind(x,y,1)

res <- construct_3d_tube(dat, radius = 0.5, len = 100)

plot3D::surf3D(res$x_mat, res$y_mat, res$z_mat, bty = "f",
               xlim = c(0,12), ylim = c(-6,6), zlim = c(-6,6))


#
# radius = 1
# len = 20
# stopifnot(nrow(dat) > 2)
#
# n <- nrow(dat)
# circle_list <- lapply(1:nrow(dat), function(x){
#   direction <- .find_adjacent_directions(dat, x)
#   res <- .find_basis_vectors(direction)
#   basis_vec1 <- res$vec1; basis_vec2 <- res$vec2
#
#   .construct_3d_circle(dat[x,], radius, basis_vec1, basis_vec2, len = len)
# })
#
# for(i in 2:length(circle_list)){
#   circle_list[[i]] <- .find_correct_orientation(circle_list[[i-1]], circle_list[[i]])
# }
#
# # form the matrices
# x_mat <- matrix(NA, nrow = nrow(circle_list[[1]]), ncol = length(circle_list))
# y_mat <- x_mat; z_mat <- x_mat
# for(i in 1:length(circle_list)){
#   x_mat[,i] <- circle_list[[i]][,1]
#   y_mat[,i] <- circle_list[[i]][,2]
#   z_mat[,i] <- circle_list[[i]][,3]
# }

