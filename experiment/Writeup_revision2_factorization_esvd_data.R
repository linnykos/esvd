rm(list=ls())
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 5,
                  2, 2, 50, 50,
                  rep(1:4, each = 8),
                  rep(c(1, 1/400, 1/350, 1/1000), each = 8),
                  rep(c(1,2, rep(3,3), rep(4,3)), times = 4),
                  rep(c(1,1, c(25, 50, 200), c(1,2,4)), times = 4),
                  rep(c(3000, rep(100, 7)), times = 4))
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "true_scalar", "true_r", "max_iter",
                        "true_distr",
                        "modifier",
                        "fitting_distr",
                        "fitting_param",
                        "max_val")
trials <- 5
ncores <- 20

################

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

vec <- paramMat[1,]
n_each <- vec["n_each"]
d_each <- vec["d_each"]
sigma <- vec["sigma"]
total <- vec["total"]
modifier <- vec["modifier"]

set.seed(10)
res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)

png(filename = "../figure/experiment/Revision_writeup2_simulation_data.png", height = 1000, width = 2500, res = 300,
    units = "px")
par(mfrow = c(1,3))
plot(res$cell_mat[,1], res$cell_mat[,2], asp = T, col = rep(1:4, each = n_each), pch = 16,
     xlab = "Latent dimension 1", ylab = "Latent dimension 2", main = "True cell latent positions")
legend("topright", c("Cell type 1", "Cell type 2", "Cell type 3", "Cell type 4"),
       bty = "n", fill = 1:4)

plot(res$gene_mat[,1], res$gene_mat[,2], asp = T, pch = 16,
     xlab = "Latent dimension 1", ylab = "Latent dimension 2", main = "True gene latent positions")

prod_mat <- res$cell_mat %*% t(res$gene_mat)
# clockwise90 = function(a) { t(a[nrow(a):1,]) } # Handy rotate function
# image(clockwise90(prod_mat), asp = 200/240)
hist(as.numeric(prod_mat), breaks = 20, col = "gray", main = "Values of natural parameter matrix",
     xlab = "Value")
graphics.off()
