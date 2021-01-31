rm(list=ls())
load("../results/step5_zeisel_comparison.RData")
library(Seurat)

cluster_labels <- as.numeric(label_vec);
levels(label_vec)
color_func <- function(alpha = 0.2){
  c(grDevices::rgb(240/255, 228/255, 66/255, alpha), #yellow
    grDevices::rgb(86/255, 180/255, 233/255, alpha), #skyblue
    grDevices::rgb(0/255, 158/255, 115/255, alpha), #bluish green
    grDevices::rgb(0/255, 114/255, 178/255,alpha), #blue
    grDevices::rgb(230/255, 159/255, 0/255,alpha), #orange
    grDevices::rgb(100/255, 100/255, 100/255, alpha), #gray
    grDevices::rgb(238/255, 31/255, 239/255, alpha))  #magenta
}
col_vec <- color_func(1)

method_vec <- c("eSVD", "SVD", "ICA", "Diffusion map", "NMF", "Isomap", "ZINB-WaVE",
                "pCMF")
fit_list <- list(esvd_fit, svd_embedding, ica_fit, diffusion_fit,
                 nnmf_fit, isomap_fit, zinbwave_fit, pcmf_fit)

umap_list <- vector("list", 8)
for(i in 1:length(fit_list)){
  set.seed(10); umap_list[[i]] <- Seurat::RunUMAP(fit_list[[i]])@cell.embeddings
}

order_vec <- c(1,7,8,2,5,3,6,4)
umap_list <- umap_list[order_vec]; method_vec <- method_vec[order_vec]

# try it our embedding
grDevices::png(filename = paste0("../../esvd_results/figure/main/zeisel_comparison.png"),
               height = 2200, width = 2000, res = 300,
               units = "px")
par(mfrow = c(3,3), mar = c(4,4,2,0.5))
for(i in 1:8){
  graphics::plot(umap_list[[i]][,1], umap_list[[i]][,2],
       col = col_vec[cluster_labels], pch = 16, asp = T,
       main = "", xlab = "UMAP dim. 1", ylab = "UMAP dim. 2")
  graphics::title(method_vec[i], line = 0.3)
}
graphics.off()

################################################

purity_list <- vector("list", 8)
neigh_size_vec <- rep(NA, 8)
for(i in 1:8){
  print(i)
  neigh_size_vec[i] <- determine_minimium_neighborhood_size(fit_list[[i]])
  set.seed(10)
  purity_list[[i]] <- unlist(compute_purity(fit_list[[i]], cluster_labels,
                        neighborhood_size = neigh_size_vec[i], num_samples = 5000)$value_list)
}

sapply(purity_list[order_vec], function(x){quantile(x, probs = c(0.05,0.1,0.15))})

