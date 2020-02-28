rm(list=ls())
library(profvis)

# set.seed(10)
# cell_pop <- matrix(c(4,10, 25,100,
#                      60,80, 25,100,
#                      40,10, 60,80,
#                      60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
# h <- nrow(cell_pop)
# n_each <- 1000
# dat <- do.call(rbind, lapply(1:h, function(x){
#   pos <- stats::runif(n_each)
#   cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.5),
#         pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.5))
# }))
# cluster_labels <- rep(1:4, each = n_each)

set.seed(10)
cluster_labels <- sample(1:4, 1000, replace = T)
dat <- MASS::mvrnorm(1000, rep(0, 5), diag(5))

plot(dat[,1], dat[,2], asp = T, col = cluster_labels, pch = 16)

p <- profvis::profvis({
  slingshot(dat, cluster_labels, starting_cluster = 1)
})

p

#############

# mat <- matrix(1:30,6,5)
# weight <- c(1:6)
# weight*mat
# diag(weight)%*%mat
