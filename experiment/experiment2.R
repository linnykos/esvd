rm(list=ls())
library(eSVD)
source("../simulation/factorization_generator.R")

paramMat <- cbind(c(10, 50, 100), 120, 5,
                  2, 50, 1/250, 1000,
                  50)
colnames(paramMat) <- c("n_each", "d_each", "sigma",
                        "k", "max_iter", "modifier", "max_val",
                        "size")

col_func2 <- function(alpha){
  c( rgb(86/255, 180/255, 233/255, alpha), #skyblue
     rgb(240/255, 228/255, 66/255, alpha), #yellow
     rgb(0/255, 158/255, 115/255, alpha), #bluish green
     rgb(230/255, 159/255, 0/255,alpha)) #orange
}
col_vec <- col_func2(1)

cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)

# compute the grid of the density
set.seed(10)
n_pop <- 500
vec <- paramMat[1,]
pop_res <- generate_natural_mat(cell_pop, gene_pop, n_pop, vec["d_each"], 0.01, vec["modifier"])

plot(pop_res$cell_mat[,1], pop_res$cell_mat[,2], asp = T, pch = 16, col = rep(1:4, each = n_pop))

#####

dat = pop_res$cell_mat
cluster_labels = rep(1:4, each = n_pop)
starting_cluster = 1
cluster_group_list = NA
use_initialization = F

stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])
stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0), length(unique(cluster_labels)) == max(cluster_labels))
if(all(!is.na(cluster_group_list))){
  tmp <- unlist(cluster_group_list)
  stopifnot(length(tmp) == length(unique(tmp)), length(tmp) == length(unique(cluster_labels)))
}

### construct the distance matrix
dist_mat <- .compute_cluster_distances(dat, cluster_labels)
dist_mat <- dist_mat^2

###############

cluster_labels = rep(1:4, each = n_pop)
par(mfrow = c(1,4))
for(i in 1:4){
  idx <- which(cluster_labels == i)
  plot(dat[idx,1], dat[idx,2], xlim = range(dat[,1]), ylim = range(dat[,2]), col = i, pch = 16, asp = T)
}

.compute_cluster_distances(dat, cluster_labels)
