rm(list=ls())
load("../results/step4_factorization.RData")
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))
alpha_val <- 0.3
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]

num_cell <- length(unique(cell_type_vec))

combn_mat <- combn(4,2)
for(x in 1:ncol(combn_mat)){
  i1 <- combn_mat[1,x]; i2 <- combn_mat[2,x]
  xlim <- range(res_our$u_mat[,i1])
  ylim <- range(res_our$u_mat[,i2])
  order_vec <- c(6,5,1,4,2,3) #these are in alphabetical order from "CO", "MF", "MO", "NF", "OP", "PP"
  name_vec <- c("Pdgfra+ (1)", "Precusor (2)", "Differentiated-commited\nprecusor (3)", "Newly formed (4)",
                "Myelin-forming (5)", "Mature (6)")

  png(paste0("../figure/experiment/Week35_marques_latent_", i1, i2, ".png"), height = 1500, width = 2000, res = 300, units = "px")
  par(mfrow = c(2,3), mar = c(4,4,4,0.5))
  for(i in 1:6){
    idx <- which(as.numeric(as.factor(cell_type_vec)) == order_vec[i])
    plot(res_our$u_mat[-idx,i1], res_our$u_mat[-idx,i2], asp = T, pch = 16,
         col = rgb(0.8, 0.8, 0.8), xlim = xlim, ylim = ylim,
         main = name_vec[i], xlab = paste0("Latent dimension ", i1),
         ylab = paste0("Latent dimension ", i2))

    lines(c(-100,100), rep(0,2), lwd = 2, lty = 2)
    lines(rep(0,2), c(-100,100), lwd = 2, lty = 2)

    points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
           col = rgb(1,1,1), cex = 1.5)

    points(res_our$u_mat[idx,i1], res_our$u_mat[idx,i2], pch = 16,
           col = col_idx[idx], cex = 1.5)
  }
  graphics.off()
}

#############

# let's try a 3d plot
car::scatter3d(x = res_our$u_mat[,1], y = res_our$u_mat[,2], z = res_our$u_mat[,3],
               surface=FALSE, groups = as.factor(cell_type_vec))

