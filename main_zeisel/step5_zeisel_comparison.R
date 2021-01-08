rm(list=ls())
load("../results/step3_zeisel_factorization.RData")
session_info <- sessionInfo(); date_of_run <- Sys.time()

library(NMF)
library(dimRed)
library(destiny)
library(umap)
library(pCMF)
library(SummarizedExperiment)
library(zinbwave)
library(fastICA)

################################################

# ICA
print("Starting ICA")
set.seed(10)
dimRed_obj <- dimRed::dimRedData(dat)
fastica_obj <- dimRed::FastICA()
fastica_obj@stdpars$ndim <- 20
suppressMessages(emb <- fastica_obj@fun(dimRed_obj, pars = fastica_obj@stdpars))
ica_fit <- emb@data@data

save.image("../results/step5_zeisel_comparison.RData")

# diffusion
print("Starting diffusion map")
set.seed(10)
tmp <- destiny::DiffusionMap(dat, n_pcs = 20)
diffusion_fit <- tmp@eigenvectors[,1:20] %*% diag(tmp@eigenvalues[1:20])

save.image("../results/step5_zeisel_comparison.RData")

# nnmf
print("Starting NNMF")
set.seed(10)
dimRed_obj <- dimRed::dimRedData(dat)
nnmf_obj <- dimRed::NNMF()
nnmf_obj@stdpars$ndim <- 20
suppressMessages(emb <- nnmf_obj@fun(dimRed_obj, pars = nnmf_obj@stdpars)) # takes around 30 min
nnmf_fit <- emb@data@data

save.image("../results/step5_zeisel_comparison.RData")

# isomap
print("Starting Isomap")
set.seed(10)
dimRed_obj <- dimRed::dimRedData(dat)
isomap_obj <- dimRed::Isomap()
isomap_obj@stdpars$ndim <- 20
isomap_obj@stdpars$knn <- round(nrow(dat)/10)
suppressMessages(emb <- isomap_obj@fun(dimRed_obj, pars = isomap_obj@stdpars))
isomap_fit <- emb@data@data

save.image("../results/step5_zeisel_comparison.RData")

# zinbwave
print("Starting Zinbwave")
load("../results/step0_zeisel_preprocessing.RData")
var_idx <- idx
load("../results/step3_zeisel_factorization.RData")

set.seed(10)
dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_count[,var_idx])))
tmp <- zinbwave::zinbwave(dat_se, K = 20, maxiter.optimize = 100, normalizedValues = F,
                          verbose = T,
                          commondispersion = F, BPPARAM = BiocParallel::SerialParam())
zinbwave_fit <- SingleCellExperiment::reducedDims(tmp)$zinbwave

save.image("../results/step5_zeisel_comparison.RData")

# pcmf
print("Starting PCMF")
set.seed(10)
tmp <- pCMF::pCMF(dat_count[,var_idx], K = 20, sparsity = F, verbose = T, ncores = 1)
pcmf_fit <- tmp$factor$U

esvd_fit <- esvd_embedding$u_mat

keep_vars <- c("ica_fit", "diffusion_fit", "nnmf_fit", "isomap_fit", "zinbwave_fit",
               "pcmf_fit", "esvd_fit", "label_vec", "svd_embedding", "dat",
               "session_info", "date_of_run")
var_list <- ls()
var_list <- var_list[!var_list %in% keep_vars]
rm(list = var_list)

source_code_info <- readLines("../main_zeisel/step5_zeisel_comparison.R")
save.image("../results/step5_zeisel_comparison.RData")












