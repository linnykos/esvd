rm(list=ls())
load("../results/step0_screening.RData")

res_hvg <- descend::findHVG(res_descend, threshold = 50)
length(res_hvg$HVG.genes)

spca_summary <- cbind(v_seq, t(sapply(res_list, function(x){c(length(unique(sort(unlist(apply(x$v, 2, function(y){which(y != 0)}))))), x$prop.var.explained[5])})))
idx <- min(which(spca_summary[,2] == max(spca_summary[,2])))
target_var <- spca_summary[idx,3]
idx <- min(intersect(which(spca_summary[,2] >= 500), which(spca_summary[,3] >= 0.9*target_var)))
spca_idx <- sort(unlist(apply(res_list[[idx]]$v, 2, function(x){which(x != 0)})))

descend_idx <- which(colnames(dat) %in% res_hvg$HVG.genes)
gene_idx <- sort(unique(c(spca_idx, descend_idx)))
dat <- dat[,gene_idx]

dat_impute <- dat

###########


dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat_impute)))
tmp <- zinbwave::zinbwave(dat_se, K = 5, maxiter.optimize = 100, normalizedValues = T,
                          commondispersion = F)
save.image("../experiment/Writeup_revision7_zinbwave.RData")

#######################

zinbwave_embedding <- SingleCellExperiment::reducedDims(tmp)$zinbwave
color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})


num_order_vec <- c(5, rep(3,2), c(6,1,1,4,6,4), rep(2,2),  rep(5,2))
col_vec <- color_func(1)[num_order_vec]
col_vec3 <- color_func(0.5)[num_order_vec]
col_name <- c("orange", rep("bluish green", 2), c("bluish green", "yellow", "yellow", "blue", "bluish green", "blue"), rep("skyblue", 2), rep("orange", 2))
order_vec <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       order = order_vec,
                       col_name = col_name,
                       col_code = col_vec)
col_info

png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup7_main_zinbwave.png"),
    height = 1250, width = 3500, res = 300,
    units = "px")
par(mar = c(4,4,4,1), mfrow = c(1,3))
plot(zinbwave_embedding[,1], zinbwave_embedding[,2],
     col = col_vec[cluster_labels], pch = 16, asp = T,
     main = "ZINB-WaVE embedding", xlab = "Latent dimension 1", ylab = "Latent dimension 2")
plot(zinbwave_embedding[,1], zinbwave_embedding[,3],
     col = col_vec[cluster_labels], pch = 16, asp = T,
     main = "ZINB-WaVE embedding", xlab = "Latent dimension 1", ylab = "Latent dimension 3")
plot(zinbwave_embedding[,2], zinbwave_embedding[,3],
     col = col_vec[cluster_labels], pch = 16, asp = T,
     main = "ZINB-WaVE embedding", xlab = "Latent dimension 2", ylab = "Latent dimension 3")
graphics.off()
