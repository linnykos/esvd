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

# let's try to make a tube?
seq_vec <- seq(0, 2*pi, length.out = 20)
x <- cos(seq_vec); y <- sin(seq_vec)
x_mat <- do.call(cbind, lapply(1:10, function(i){x}))
y_mat <- do.call(cbind, lapply(1:10, function(i){y}))
z_mat <- do.call(rbind, lapply(1:20, function(i){1:10}))

surf3D(x_mat, y_mat, z_mat, colvar = z_mat, bty = "f")
