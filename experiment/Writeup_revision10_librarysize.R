rm(list=ls())
load("../../raw_data/marques.RData")

dat <- marques$counts
colnames(dat) <- gsub("_", "-", colnames(dat))
dat_count <- dat

cell_type_vec <- as.character(marques$cell.info$cell.type)
cell_type_vec <- as.factor(cell_type_vec)
dim(dat)

library_size <- rowSums(dat)

cluster_labels <- as.numeric(cell_type_vec)
custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))

tmp_df <- data.frame(val = library_size,  type = sapply(1:13, function(y){which(sapply(custom_cluster_group_list, function(x){y %in% x}))})[cluster_labels])

vioplot::vioplot(tmp_df$val[tmp_df$type == 1],
                 tmp_df$val[tmp_df$type == 2],
                 tmp_df$val[tmp_df$type == 3],
                 tmp_df$val[tmp_df$type == 4],
                 tmp_df$val[tmp_df$type == 5],
                 tmp_df$val[tmp_df$type == 6])
