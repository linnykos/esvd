rm(list=ls())
max_val = 50
family = "neg_binom"
width = 0.8
scalar = 10
angle_val = 45

seq_vec <- seq(0, max_val, length.out = 500)

interval_mat <- sapply(seq_vec, function(x){
  .compute_prediction_interval_from_mean(x, family = family, width = width, scalar = scalar)
})

principal_line <- seq_vec * tan(angle_val*pi/180)

bool_vec <- apply(cbind(interval_mat[1,] <= principal_line, interval_mat[2,] >= principal_line), 1, all)
bool <- sum(bool_vec)/length(bool_vec) >= 0.95

plot(seq_vec, interval_mat[2,])
points(seq_vec, interval_mat[1,], col = "red")
points(seq_vec, principal_line, col = "blue")
