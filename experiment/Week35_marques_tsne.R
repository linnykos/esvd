rm(list=ls())
load("../results/step4_factorization_spca.RData")

set.seed(10)
tsne_res <- tsne::tsne(dat_impute, perplexity = 50, k = 2)

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
# keep only the first two letters
cell_type_vec <- as.character(sapply(cell_type_vec, function(x){substr(x,1,2)}))

alpha_val <- 1
col_idx <- c(rgb(238/255,204/255,17/255,alpha_val), #goldenrod
             rgb(227/255,73/255,86/255,alpha_val), #red
             rgb(100/255,140/255,252/255,alpha_val), #blue
             rgb(129/255,199/255,124/255,alpha_val), #green
             rgb(238/255,204/255,17/255,alpha_val),
             rgb(238/255,204/255,17/255,alpha_val))[as.numeric(as.factor(cell_type_vec))]


png("../figure/experiment/Week35_marques_tsne.png", height = 1200, width = 1200, res = 300, units = "px")
plot(tsne_res[,1], tsne_res[,2], pch = 16, col = col_idx, asp = T,
     xlab = "Latent dimension 1", ylab = "Latent dimension 2")
graphics.off()
