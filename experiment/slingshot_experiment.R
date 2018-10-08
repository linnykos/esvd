rm(list=ls())
source("../experiment/Week23_simulation_generator.R")

set.seed(10)
simulation <- .data_generator(total = 200)
dat <- simulation$cell_mat
cluster_labels <- rep(1:simulation$h, each = simulation$n_each)

omega <- Inf
slingshot_lineage <- slingshot::getLineages(data = simulation$cell_mat_org,
                                      clusterLabels = rep(1:simulation$h, each = simulation$n_each),
                                      omega = omega)
slingshot_res <- slingshot::slingshot(data = simulation$cell_mat_org,
                                            clusterLabels = rep(1:simulation$h, each = simulation$n_each),
                                            omega = omega)

slingshot_lineage@adjacency

#######

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green
plot(simulation$cell_mat_org[,1], simulation$cell_mat_org[,2],
     xlim = range(c(simulation$cell_mat_org[,1], 0)),
     ylim = range(c(simulation$cell_mat_org[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(slingshot_res)

###########

res <- .get_lineages(dat, cluster_labels)
lineages <- res$lineages
cluster_mat <- res$cluster_mat

shrink = TRUE; extend = 'y'; reweight = TRUE; reassign = TRUE; thresh = 0.001
maxit = 15; stretch = 2; smoother = 'smooth.spline'; shrink_method = 'cosine'
allow.breaks = TRUE
