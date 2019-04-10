rm(list=ls())
load("../results/step5_clustering_spca.RData")

# let's start with some basic investigation
our_sd_val$sd_val
naive_sd_val$sd_val

naive_sd <- max(sapply(naive_sd_val$mat_list, function(x){
  quantile(apply(x, 2, function(x){
    quantile(x,probs = 0.95)
  }), probs = 0.95)
}))

our_sd <- max(sapply(our_sd_val$mat_list, function(x){
  quantile(apply(x, 2, function(x){
    quantile(x,probs = 0.95)
  }), probs = 0.95)
}))

###################################

# basic 3d plots
# colors from http://mkweb.bcgsc.ca/colorblind/img/colorblindness.palettes.trivial.png
col_func <- function(alpha){
  col_vec <- numeric(length(unique(cluster_labels)))
  for(i in 1:3){
    col_vec[cluster_group_list[[i]]] <- rgb(240/255, 228/255, 66/255,alpha) #yellow
  }
  col_vec[cluster_group_list[[4]]] <- rgb(0/255, 158/255, 115/255,alpha) #bluish green
  col_vec[cluster_group_list[[5]]] <- rgb(230/255, 159/255, 0/255,alpha) #red
  col_vec[cluster_group_list[[6]]] <- rgb(86/255, 180/255, 233/255,alpha) #blue

  col_vec[c(1,2,4,10,12)]
}

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

slingshot_3dplot(res_our$u_mat[,1:3], cluster_labels,
                 bg_col_vec = col_func(0.1), cluster_col_vec = col_func(1),
                 breaks = c(0, 1.5, 3.5, 9.5, 11.5, 14), curves = our_curves$curves,
                 pch = 16, lwd = 2, main = "Our lineages",
                 xlab = "Latent dimension 1",
                 ylab = "Latent dimension 2",
                 zlab = "Latent dimension 3")

# try a grid of angles
paramMat <- as.matrix(expand.grid(seq(0, 360, length.out = 9),
                                  seq(0, 360, length.out = 9)))
sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')

  png(paste0("../figure/experiment/Writeup22_3dplots/Writeup22_our_lineage_theta",
             paramMat[x,1], "_phi", paramMat[x,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0.5, 0.5, 4, 0.5))
  slingshot_3dplot(res_our$u_mat[,1:3], cluster_labels,
                   bg_col_vec = col_func(0.1), cluster_col_vec = col_func(1),
                   breaks = c(0, 1.5, 3.5, 9.5, 11.5, 14), curves = our_curves$curves,
                   pch = 16, lwd = 2, main = "Our lineages",
                   theta = paramMat[x,1], phi = paramMat[x,2],
                   xlab = "Latent dimension 1",
                   ylab = "Latent dimension 2",
                   zlab = "Latent dimension 3")
  graphics.off()
})


