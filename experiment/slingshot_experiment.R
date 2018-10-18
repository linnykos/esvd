rm(list=ls())
source("../experiment/Week25_simulation_generator.R")

set.seed(10)
simulation <- .data_generator(total = 200)
dat <- simulation$cell_mat
cluster_labels <- rep(1:simulation$h, each = simulation$n_each)

res <- .get_lineages(dat, cluster_labels)
lineages <- res$lineages
lineages[[1]] <- c("1","3","2")
lineages[[2]] <- c("1","3","4")
cluster_mat <- res$cluster_mat
curves <- .get_curves(dat, cluster_mat, lineages, reassign = F)

#######

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green
plot(simulation$cell_mat[,1], simulation$cell_mat[,2],
     xlim = range(c(simulation$cell_mat[,1], 0)),
     ylim = range(c(simulation$cell_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")

curves <- pcurves
for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  lines(curves[[i]]$s[ord,1], curves[[i]]$s[ord,2])
}

# lines(s[,1], s[,2])

############################


shrink = TRUE; extend = 'y'; reweight = TRUE; reassign = F; thresh = 0.001
maxit = 15; stretch = 2; smoother = 'smooth.spline'; shrink_method = 'cosine'
allow.breaks = TRUE

#####

shrink <- as.numeric(shrink)
# CHECKS
if(shrink < 0 | shrink > 1){
  stop("'shrink' parameter must be logical or numeric between",
       "0 and 1")
}

smoother_func <- .smoother_slingshot(smoother)

### setup
num_lineage <- length(grep("Lineage", names(lineages))) # number of lineages
clusters <- colnames(cluster_mat)
n <- nrow(dat)
p <- ncol(dat)
nclus <- length(clusters)
centers <- .compute_clustercenter(dat, cluster_mat)

if(p == 1){
  centers <- t(centers)
  rownames(centers) <- clusters
}

rownames(centers) <- clusters
W <- vapply(seq_len(num_lineage), function(l){
  rowSums(cluster_mat[, lineages[[num_lineage]], drop = FALSE])
}, rep(0,nrow(dat))) # weighting matrix
rownames(W) <- rownames(dat)
colnames(W) <- names(lineages)[seq_len(num_lineage)]
W.orig <- W
D <- W; D[,] <- NA

### determine curve hierarchy
C <- as.matrix(vapply(lineages[seq_len(num_lineage)], function(lin) {
  vapply(clusters, function(cluster_id) {
    as.numeric(cluster_id %in% lin)
  }, 0)
}, rep(0,nclus)))
rownames(C) <- clusters
segmnts <- unique(C[rowSums(C)>1,,drop = FALSE])
segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE), ,
                   drop = FALSE]
avg_order <- list() #avg_order will have length equal to number of connected components, i think
for(i in seq_len(nrow(segmnts))){
  idx <- segmnts[i,] == 1
  avg_order[[i]] <- colnames(segmnts)[idx]
  new_col <- rowMeans(segmnts[,idx, drop = FALSE])
  segmnts <- cbind(segmnts[, !idx, drop = FALSE], new_col)
  colnames(segmnts)[ncol(segmnts)] <- paste('average',i,sep='')
}

### initial curves are piecewise linear paths through the tree
pcurves <- list()
for(l in seq_len(num_lineage)){
  idx <- which(W[,l] > 0)
  line_initial <- centers[clusters %in% lineages[[l]], ,
                          drop = FALSE]
  line_initial <- line_initial[match(lineages[[l]],
                                     rownames(line_initial)),  ,
                               drop = FALSE]
  K <- nrow(line_initial)

  curve <- princurve::project_to_curve(dat[idx, ,drop = FALSE],
                                       s = line_initial, stretch = 9999) #note: this changes line_initial
  curve$dist_ind <- abs(curve$dist_ind)

  sample_idx <- sort(unique(unlist(lapply(as.numeric(lineages[[l]]), function(x){which(cluster_mat[,x] == 1)}))))
  pcurve <- princurve::project_to_curve(dat[sample_idx,], s = curve$s[curve$ord, ,drop=FALSE],
                                        stretch=0)
  # force non-negative distances
  pcurve$dist_ind <- abs(pcurve$dist_ind)
  # force pseudotime to start at 0
  pcurve$lambda <- pcurve$lambda - min(pcurve$lambda,
                                       na.rm=TRUE)
  pcurve$w <- W[sample_idx,l]
  pcurves[[l]] <- pcurve
  D[sample_idx,l] <- abs(pcurve$dist_ind)
}

### track distances between curves and data points to determine convergence
dist_new <- sum(abs(D[W>0]), na.rm=TRUE)

it <- 0
hasConverged <- FALSE

it <- it + 1
dist_old <- dist_new

if(reweight | reassign){
  ordD <- order(D)
  W_prob <- W/rowSums(W)
  WrnkD <- cumsum(W_prob[ordD]) / sum(W_prob) #ERROR?
  Z <- D
  Z[ordD] <- WrnkD #ERROR?
}

if(reweight){
  Z_prime <- 1-Z^2
  Z_prime[W==0] <- NA
  W0 <- W
  W <- Z_prime / matrixStats::rowMaxs(Z_prime, na.rm = TRUE) #rowMins(D) / D
  W[is.nan(W)] <- 1 # handle 0/0
  W[is.na(W)] <- 0
  W[W > 1] <- 1
  W[W < 0] <- 0
  W[W0==0] <- 0
}

if(reassign){
  # add if z < .5
  idx <- Z < .5
  W[idx] <- 1 #(rowMins(D) / D)[idx]

  # drop if z > .9 and w < .1
  ridx <- matrixStats::rowMaxs(Z, na.rm = TRUE) > .9 &
    matrixStats::rowMins(W, na.rm = TRUE) < .1
  W0 <- W[ridx, ]
  Z0 <- Z[ridx, ]
  W0[!is.na(Z0) & Z0 > .9 & W0 < .1] <- 0
  W[ridx, ] <- W0
}

for(l in seq_len(num_lineage)){
  pcurve_safe <- pcurves[[l]]
  pcurve <- pcurves[[l]]
  s <- pcurve$s
  ordL <- order(pcurve$lambda) #Q: why are there lambdas even for those not on the curve?
  sample_idx <- sort(unique(unlist(lapply(as.numeric(lineages[[l]]), function(x){which(cluster_mat[,x] == 1)}))))
  s <-  .smoother_func_better(pcurve$lambda, dat[sample_idx,])
  new_pcurve <- princurve::project_to_curve(dat[sample_idx,], s = s, stretch = stretch)
  new_pcurve$dist_ind <- abs(new_pcurve$dist_ind)
  new_pcurve$lambda <- new_pcurve$lambda -
    min(new_pcurve$lambda, na.rm = TRUE)
  new_pcurve$w <- W[,l]
  pcurves[[l]] <- new_pcurve
}

