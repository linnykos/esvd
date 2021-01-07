rm(list=ls())
load("../results/step3_zeisel_factorization.RData")
cluster_labels <- as.numeric(label_vec)

library(NMF)

################################################

neigh_size <- determine_minimium_neighborhood_size(esvd_embedding$u_mat)
set.seed(10)
res <- compute_purity(esvd_embedding$u_mat, cluster_labels,
               neighborhood_size = neigh_size, num_samples = 5000)

vec <- unlist(res_list[[1]]$value_list)
mean(vec)
mean(sort(vec)[1:round(length(vec)/4)])
quantile(sort(vec)[1:round(length(vec)/4)], probs = c(0, 0.25, 0.5, 0.75, 1))

# try it our embedding
set.seed(10)
zz <- Seurat::RunUMAP(esvd_embedding$u_mat)@cell.embeddings
plot(zz[,1], zz[,2], col = cluster_labels, pch = 16, asp = T)

########

# SVD
zz <- RSpectra::svds(dat, k = 20)
tmp <- zz$u %*% diag(sqrt(zz$d))
set.seed(10)
zz <- Seurat::RunUMAP(tmp)@cell.embeddings
plot(zz[,1], zz[,2], col = cluster_labels, pch = 16, asp = T)

#####

# ICA
set.seed(10)
dimRed_obj <- dimRed::dimRedData(dat)
fastica_obj <- dimRed::FastICA()
fastica_obj@stdpars$ndim <- 20
suppressMessages(emb <- fastica_obj@fun(dimRed_obj, pars = fastica_obj@stdpars))
fit <- emb@data@data

set.seed(10)
zz <- Seurat::RunUMAP(fit)@cell.embeddings
plot(zz[,1], zz[,2], col = cluster_labels, pch = 16, asp = T)

neigh_size <- determine_minimium_neighborhood_size(fit)
set.seed(10)
res <- compute_purity(fit, cluster_labels,
                      neighborhood_size = neigh_size, num_samples = 5000)
vec <- res$value_list[[1]]
mean(vec)
mean(sort(vec)[1:round(length(vec)/4)])
quantile(sort(vec)[1:round(length(vec)/4)], probs = c(0, 0.25, 0.5, 0.75, 1))

####

# diffusion
set.seed(10)
tmp <- destiny::DiffusionMap(dat, n_pcs = 20)
fit <- tmp@eigenvectors[,1:20] %*% diag(tmp@eigenvalues[1:20])

set.seed(10)
zz <- Seurat::RunUMAP(fit)@cell.embeddings
plot(zz[,1], zz[,2], col = cluster_labels, pch = 16, asp = T)

neigh_size <- determine_minimium_neighborhood_size(fit)
set.seed(10)
res <- compute_purity(fit, cluster_labels,
                      neighborhood_size = neigh_size, num_samples = 5000)
vec <- res$value_list[[1]]
mean(vec)
mean(sort(vec)[1:round(length(vec)/4)])
quantile(sort(vec)[1:round(length(vec)/4)], probs = c(0, 0.25, 0.5, 0.75, 1))

########

set.seed(10)
dimRed_obj <- dimRed::dimRedData(dat)
nnmf_obj <- dimRed::NNMF()
nnmf_obj@stdpars$ndim <- 20
suppressMessages(emb <- nnmf_obj@fun(dimRed_obj, pars = nnmf_obj@stdpars))
fit <- emb@data@data

