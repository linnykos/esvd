rm(list=ls())
dat <- read.table("../data-raw/Zhang_marker_genes.txt", sep = ",", strip.white = T)
# from figure 4 of https://www.jneurosci.org/content/jneuro/34/36/11929.full.pdf

new_mat <- do.call(rbind, lapply(1:nrow(dat), function(i){
  vec <- as.character(as.factor(dat[i,]))
  cbind(rep(vec[1], length(vec)-1), vec[-1])
}))

zhang_genes <- new_mat
zhang_genes <- as.data.frame(zhang_genes)
colnames(zhang_genes) <- c("cell_type", "enriched_genes")

usethis::use_data(zhang_genes, overwrite = T)


