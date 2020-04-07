rm(list=ls())
library(simulation)
library(eSVD)
source("../simulation/factorization_generator.R")
source("../simulation/factorization_methods.R")

paramMat <- cbind(c(10, 50, 100), 120, 5,
                  2, 50, 1/250, 1000,
                  50)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "max_iter", "modifier", "max_val",
                        "size")

trials <- 1
ncores <- NA

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

  list(fit = fit)
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

.compute_true_density <- function(cell_mat, grid_size,
                                  sigma = 0.05){
  spacing <- 0.25
  xrange <- c(floor(min(cell_mat[,1]-.5)/spacing)*spacing, ceiling(max(cell_mat[,1]+.5)/spacing)*spacing)
  yrange <- c(floor(min(cell_mat[,2]-.5)/spacing)*spacing, ceiling(max(cell_mat[,2]+.5)/spacing)*spacing)

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
pop_res <- generate_natural_mat(cell_pop, gene_pop, n_pop, 500, 0)
den_res <- .compute_true_density(pop_res$cell_mat, 151, 0.1)
mat <- den_res$density_mat
rownames(mat) <- as.numeric(rownames(mat))
colnames(mat) <- as.numeric(colnames(mat))

# compute the points associated with each cluster
## identify all the high-probability regions
idx <- which(mat > quantile(mat, probs = 0.95), arr.ind = T)
## translate into coordinates
idx[,1] <- as.numeric(rownames(mat))[idx[,1]]
idx[,2] <- as.numeric(colnames(mat))[idx[,2]]
idx <- idx[,c(2,1)]
## assign each point to clusters
cluster_identifier <- as.numeric(apply(idx, 1, function(x){
  i <- which.min(abs(apply(pop_res$cell_mat, 1, function(y){.l2norm(x-y)})))
  floor((i-1)/n_pop)+1
}))

# compute the population lineage
pop_lineage <- eSVD::slingshot(pop_res$cell_mat,
                                     rep(1:4, each = n_pop),
                                     starting_cluster = 1, verbose = F, reduction_percentage = 0.1)

png("../figure/simulation/example_trajectories.png", height = 960, width = 2500, res = 300, units = "px")
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

axis(1, at = seq(-3,0,by=1), labels = T, las=2)
axis(2, at = seq(-3,2,by=1), labels = T, las=2)

for(j in 1:length(pop_lineage$curves)){
  ord <- pop_lineage$curves[[j]]$ord
  lines(pop_lineage$curves[[j]]$s[ord,1],
        pop_lineage$curves[[j]]$s[ord,2], lwd = 2)
}

# plot the centers of each cluster
pop_means <- t(sapply(1:4, function(i){
  colMeans(idx[which(cluster_identifier == i),])
}))
points(pop_means[,1], pop_means[,2], pch = 21, cex = 2, bg = col_vec)

# HOT FIX
for(i in 1:length(res)){
  tmp <- res[[i]][[1]]$res_our; tmp[,1] <- -tmp[,1]
  n <- paramMat[i,"n_each"]
  samp_cluster_identifier <- rep(1:4, each = n)

  # compute the rescaling factor
  samp_means <- t(sapply(1:4, function(i){
    colMeans(tmp[which(samp_cluster_identifier == i),])
  }))
  rescale_factor <- median(pop_means/samp_means)
  tmp <- tmp*rescale_factor

  # sig <- ifelse(mean(tmp[1:n_seq[i],2]) > mean(tmp[(3*n_seq[i]+1):(4*n_seq[i]),2]), -1, 1)
  # tmp[,2] <- sig*tmp[,2]
  plot(tmp[,1], tmp[,2],
       pch = 16, col = col_vec[rep(1:4, each = n)],
       asp = T,
       xlim = range(as.numeric(colnames(mat))),
       ylim = range(as.numeric(rownames(mat))),
       xlab = "Latent dimension 1", ylab = "Latent dimension 2", axes = F, cex.lab = 1.25,
       main = paste0("Estimated embedding\n(n = ", 4*n, ")"))

  axis(1, at = seq(-3,0,by=1), labels = T, las=2)
  axis(2, at = seq(-3,2,by=1), labels = T, las=2)

  # for(j in 1:length(res[[i]][[1]]$curves_our$curves)){
  #   ord <- res[[i]][[1]]$curves_our$curves[[j]]$ord
  #   lines(-res[[i]][[1]]$curves_our$curves[[j]]$s[ord,1]*rescale_factor,
  #         res[[i]][[1]]$curves_our$curves[[j]]$s[ord,2]*rescale_factor, lwd = 2)
  # }
}
graphics.off()
