rm(list=ls())
load("../results/step4_factorization_spca.RData")

d <- 3
dat <- res_our$u_mat[,1:d]

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
alpha_val <- 0.2
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]

num_cell <- length(unique(cell_type_vec))

combn_mat <- combn(3,2)
for(x in 1:ncol(combn_mat)){
  i1 <- combn_mat[1,x]; i2 <- combn_mat[2,x]
  xlim <- range(dat[,i1])
  ylim <- range(dat[,i2])
  order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
  name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
                "Myelin-forming (5)", "Mature (6)")

  png(paste0("../figure/experiment/Writeup22_marques_latent_our_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
  par(mfrow = c(2,3), mar = c(4,4,4,0.5))
  for(i in 1:6){
    idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
    plot(dat[-idx,i1], dat[-idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
         main = name_vec[i], xlab = paste0("Our latent dimension ", i1),
         ylab = paste0("Our latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = rgb(1,1,1), cex = 1.5)

    points(dat[idx,i1], dat[idx,i2], pch = 16,
           col = col_idx[idx], cex = 1.5)
  }
  graphics.off()

}


##################

#estimate lineage now

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

curves <- singlecell::slingshot(res_our$u_mat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                                cluster_group_list = cluster_group_list,
                                verbose = T)


######

col_func <- function(alpha){
  col_vec <- numeric(length(unique(cluster_labels)))
  for(i in 1:3){
    col_vec[cluster_group_list[[i]]] <- rgb(238/255,204/255,17/255,alpha) #goldenrod
  }
  col_vec[cluster_group_list[[4]]] <- rgb(129/255,199/255,124/255,alpha) #green
  col_vec[cluster_group_list[[5]]] <- rgb(227/255,73/255,86/255,alpha) #red
  col_vec[cluster_group_list[[6]]] <- rgb(100/255,140/255,252/255,alpha) #blue

  col_vec[c(1,2,4,10,12)]
}


plot_3d_func <- function(dat, curves, name, theta = 60, phi = 20){
  cluster_center <- .compute_cluster_center(dat, .construct_cluster_matrix(cluster_labels))

  png(paste0("../figure/experiment/Writeup22_3dplots/Writeup22_marques_", name,
             "_lineage_theta", theta, "_phi", phi, ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0.5, 0.5, 4, 0.5))
  plot3D::scatter3D(x = dat[,1], y = dat[,2], z = dat[,3],
                    surface = FALSE, colvar = cluster_labels,
                    breaks = c(0, 1.5, 3.5, 9.5, 11.5, 14),
                    col = col_func(0.1), pch = 16, cex = 0.4, colkey = F,
                    theta = theta, phi = phi, main = paste0(name, " lineages"),
                    xlab = "Latent dimension 1",
                    ylab = "Latent dimension 2",
                    zlab = "Latent dimension 3")
  plot3D::points3D(cluster_center[,1], cluster_center[,2], cluster_center[,3],
                   colvar = 1:13,
                   breaks = c(0, 1.5, 3.5, 9.5, 11.5, 14),
                   add = T, pch = 16, cex = 2, col = col_func(1), colkey = F)

  reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25
  for(k in 1:length(curves$curves)){
    ord <- curves$curves[[k]]$ord
    plot3D::lines3D(x = curves$curves[[k]]$s[ord, 1]*reduction_factor,
                    y = curves$curves[[k]]$s[ord, 2]*reduction_factor,
                    z = curves$curves[[k]]$s[ord, 3]*reduction_factor,
                    add = T, lwd = 2, colkey = F, col = "black")
  }
  graphics.off()
}

plot_3d_func(dat, curves, "Our")


paramMat <- as.matrix(expand.grid(seq(0, 360, length.out = 17),
                                  seq(0, 360, length.out = 17)))
sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')
  plot_3d_func(dat, curves, "Our", theta = paramMat[x,1], phi = paramMat[x,2])
})


