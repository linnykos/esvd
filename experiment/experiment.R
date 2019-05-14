rm(list=ls())
library(simulation)
library(singlecell)
source("../simulation/factorization_generator.R")

paramMat <- cbind(50, 120, 0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
trials <- 50
vec <- paramMat[1,]

#####################

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25)/10,
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/100, nrow = 2, ncol = 4, byrow = T)


n_each <- vec["n_each"]
d_each <- vec["d_each"]
sigma <- vec["sigma"]
total <- vec["total"]

res <- generate_natrual_mat(cell_pop, gene_pop, n_each, d_each, sigma)
nat_mat <- res$nat_mat

set.seed(10)
obs_mat <- generator_zinb_nb(nat_mat)
# quantile(obs_mat, probs = seq(0, 1, length.out=11))

# set.seed(10)
# tmp <- pCMF::pCMF(obs_mat, K = vec["k"], verbose = T, sparsity = F)
# res_pcmf <- tmp$factor$U
# plot(res_pcmf[,1], res_pcmf[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])

set.seed(10)
dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(obs_mat)))
tmp <- zinbwave::zinbwave(dat_se, K = vec["k"], maxiter.optimize = 100)
res_zinb <- tmp@reducedDims$zinbwave
plot(res_zinb[,1], res_zinb[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])


##############

set.seed(10)
init <- singlecell::initialization(obs_mat, family = "gaussian", k = vec["k"], max_val = vec["max_val"])
tmp <- singlecell::fit_factorization(obs_mat, u_mat = init$u_mat, v_mat = init$v_mat,
                                     family = "gaussian",  reparameterize = T,
                                     max_iter = 100, max_val = vec["max_val"],
                                     scalar = 1,
                                     return_path = F, cores = 1,
                                     verbose = F)
res_our <- tmp$u_mat
plot(res_our[,1], res_our[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])
cluster_labels <- rep(1:4, each = vec["n_each"])
curves_our <- singlecell::slingshot(res_our[,1:vec["k"]], cluster_labels,
                                    starting_cluster = 1,
                                    verbose = F)
curves_our$lineages
for(i in 1:3){
  ord <- curves_our$curves[[i]]$ord
  lines(curves_our$curves[[i]]$s[ord,1], curves_our$curves[[i]]$s[ord,2], lwd = 2)
}

###############

dat <- res_our[,1:vec["k"]]
starting_cluster <- 1
cluster_group_list = NA
dist_mat <- .compute_cluster_distances(dat, cluster_labels)

######
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

for(i in 1:4){
  mean_vec <- apply(dat[which(cluster_labels == i), ], 2, mean)
  cov_mat <- cov(dat[which(cluster_labels == i), ])
  epoints <- ellipse_points(mean_vec, cov_mat)
  lines(epoints[,1], epoints[,2], col = "black", lwd = 2)
}

i <- 4
dat_subset <- dat[which(cluster_labels == i), ]
mean_vec <- apply(dat_subset, 2, mean)
cov_mat <- cov(dat_subset)
epoints <- ellipse_points(mean_vec, cov_mat)
tr <- tripack::tri.mesh(epoints[-1,1], epoints[-1,2])
zz <- tripack::in.convex.hull(tr, dat_subset[,1], dat_subset[,2])
sum(zz)/length(zz)

# after intersecting, discard all other indices in the distance matrix

