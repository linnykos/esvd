rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")
source("../simulation/factorization_methods.R")

paramMat <- cbind(c(10, 50, 200), 120, 5,
                  2, 50, 1/250, 1000,
                  50)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "max_iter", "modifier", "max_val",
                        "size")

trials <- 1
ncores <- 10

################

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

rule <- function(vec){
  n_each <- vec["n_each"]
  d_each <- vec["d_each"]
  sigma <- vec["sigma"]
  modifier <- vec["modifier"]

  res <- generate_natural_mat(cell_pop, gene_pop, n_each, d_each, sigma, modifier)
  nat_mat <- res$nat_mat

  dat <- generator_esvd_nb(nat_mat, vec["size"])
  obs_mat <- round(dat * 1000/max(dat))

  list(dat = obs_mat, truth = res$cell_mat)
}

criterion <- function(dat, vec, y){
  dat_obs <- dat$dat

  set.seed(10)
  init <- eSVD::initialization(dat_obs, family = "neg_binom", k = vec["k"], max_val = 2000,
                               scalar = vec["size"])
  fit <- eSVD::fit_factorization(dat_obs, u_mat = init$u_mat, v_mat = init$v_mat,
                                 family = "neg_binom", scalar = vec["size"],
                                 max_iter = 50, max_val = 2000,
                                 return_path = F, cores = ncores, verbose = F)

  list(fit = fit, truth = dat$truth)
}

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = NA, as_list = T,
                                        filepath = "../results/illustration_example_tmp.RData",
                                        verbose = T)
save.image("../results/illustration_example.RData")

##  i <- 2; y <- 1; set.seed(y); zz1 <- criterion(rule(paramMat[i,]), paramMat[i,], y); plot(zz1$fit$u_mat[,1], zz1$fit$u_mat[,2], asp = T)

###################################
rm(list=ls())
load("../results/illustration_example.RData")
for(j in 1:3){
  print(apply(res[[j]][[1]]$fit$u_mat, 2, range))
  print(apply(res[[j]][[1]]$truth, 2, range))
}

res[[1]][[1]]$fit$u_mat <- res[[1]][[1]]$fit$u_mat %*% diag(c(-1,1))
res[[2]][[1]]$fit$u_mat <- res[[2]][[1]]$fit$u_mat %*% diag(c(1,-1))
res[[3]][[1]]$fit$u_mat <- res[[3]][[1]]$fit$u_mat %*% diag(c(-1,1))


.compute_true_density <- function(cell_mat, grid_size, xrange, yrange,
                                  sigma = 5){
  xseq <- seq(xrange[1], xrange[2], length.out = grid_size)
  yseq <- seq(yrange[1], yrange[2], length.out = grid_size)
  mat <- matrix(NA, nrow = length(yseq), ncol = length(xseq))

  colnames(mat) <- xseq; rownames(mat) <- rev(yseq)

  for(i in 1:length(xseq)){
    if(i %% floor(length(xseq)/10) == 0) cat('*')

    for(j in 1:length(yseq)){
      d <- min(apply(cell_mat, 1, function(x){
        .l2norm(x - c(xseq[i], yseq[j]))
      }))
      mat[length(yseq)-j+1,i] <- dnorm(d, sd = sigma)
    }
  }

  list(density_mat = mat)
}

.rotate <- function(a) { t(a[nrow(a):1,]) }
#############################

col_func2 <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col_vec <- col_func2(1)

# compute the grid of the density
set.seed(10)
n_pop <- 500
vec <- paramMat[1,]
pop_res <- generate_natural_mat(cell_pop, gene_pop, n_pop, vec["d_each"], 0.01, vec["modifier"])
den_res <- .compute_true_density(pop_res$cell_mat, 101, vec["sigma"]/4, xrange=  c(-1.3, 0.2), yrange = c(0.05-(2.75*1.5/2)/2, 0.05+(2.75*1.5/2)/2))
mat <- den_res$density_mat
rownames(mat) <- as.numeric(rownames(mat))
colnames(mat) <- as.numeric(colnames(mat))

# plot(pop_res$cell_mat[,1], pop_res$cell_mat[,2], asp = T, col = rep(1:4, each = n_pop), pch = 16)

# compute the population lineage
lineages <- list(Lineage1 = c(1,2,3), Lineage2 = c(1,2,4))
cluster_labels = rep(1:4, each = n_pop)
starting_cluster = 1
cluster_group_list = NA
use_initialization = F
reduction_percentage = 0.1
shrink = 1
thresh = 0.001
max_iter = 15
upscale_factor = NA
verbose = F

dat <- pop_res$cell_mat
reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*reduction_percentage
dat2 <- dat/reduction_factor

curve_res <- .get_curves(dat2, cluster_labels, cluster_group_list, lineages, shrink = shrink,
                   thresh = thresh, max_iter = max_iter, upscale_factor = upscale_factor,
                   verbose = verbose)
curves <- curve_res$pcurve_list

# adjust up
for(k in 1:length(curves)){
  curves[[k]]$s <- curves[[k]]$s*reduction_factor
}

pop_lineage <- list(lineages = lineages, curves = curves, idx = res$idx)

png("../../esvd_results/figure/experiment/example_trajectories.png", height = 960, width = 2500, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(4,4,4,0.1))

image(as.numeric(colnames(mat)),
      rev(as.numeric(rownames(mat))), .rotate(mat),
      col = grDevices::heat.colors(100, alpha = 0.5),
      xlab = "Latent dimension 1", ylab = "Latent dimension 2", asp = T,
      main = "Population trajectory", cex.lab = 1.25,
      axes = F)
contour(as.numeric(colnames(mat)),
        rev(as.numeric(rownames(mat))), .rotate(mat),
        add = T, drawlabels = F, col = rgb(0,0,0,0.5), lwd = 1,
        levels = quantile(mat, probs = c(0.25,0.5,0.75)))

axis(1)
axis(2)

for(j in 1:length(pop_lineage$curves)){
  ord <- pop_lineage$curves[[j]]$ord
  lines(pop_lineage$curves[[j]]$s[ord,1],
        pop_lineage$curves[[j]]$s[ord,2], lwd = 2)
}

# plot the centers of each cluster
pop_means <- t(sapply(1:4, function(i){
  colMeans(pop_res$cell_mat[which(cluster_labels == i),])
}))
points(pop_means[,1], pop_means[,2], pch = 21, cex = 2, bg = col_vec)

# HOT FIX
for(i in 1:length(res)){
  tmp <- res[[i]][[1]]$fit$u_mat
  n <- paramMat[i,"n_each"]
  samp_cluster_identifier <- rep(1:4, each = n)

  # compute the rescaling factor
  samp_means <- t(sapply(1:4, function(i){
    colMeans(tmp[which(samp_cluster_identifier == i),])
  }))

  # sig <- ifelse(mean(tmp[1:n_seq[i],2]) > mean(tmp[(3*n_seq[i]+1):(4*n_seq[i]),2]), -1, 1)
  # tmp[,2] <- sig*tmp[,2]
  plot(tmp[,1], tmp[,2],
       pch = 16, col = col_vec[rep(1:4, each = n)],
       asp = T,
       xlab = "Latent dimension 1", ylab = "Latent dimension 2", axes = F, cex.lab = 1.25,
       main = paste0("Estimated embedding\n(n = ", 4*n, ")"))

  axis(1)
  axis(2)

  # for(j in 1:length(res[[i]][[1]]$curves_our$curves)){
  #   ord <- res[[i]][[1]]$curves_our$curves[[j]]$ord
  #   lines(-res[[i]][[1]]$curves_our$curves[[j]]$s[ord,1]*rescale_factor,
  #         res[[i]][[1]]$curves_our$curves[[j]]$s[ord,2]*rescale_factor, lwd = 2)
  # }
}
graphics.off()
