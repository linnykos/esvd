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

# dbCon <- org.Mm.eg.db::org.Mm.eg_dbconn()
# sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
#
# vec <- sapply(1:nrow(new_mat), function(i){
#   any(c(new_mat[i,2] %in% aliasSymbol$alias_symbol, new_mat[i,2] %in% aliasSymbol$symbol))
# })
# new_mat[!vec, 2]
#
# syn_vec <- sapply(1:nrow(new_mat), function(i){
#   bool <- any(c(new_mat[i,2] %in% aliasSymbol$alias_symbol, new_mat[i,2] %in% aliasSymbol$symbol))
#   if(!bool | new_mat[i,2] %in% aliasSymbol$symbol) return(new_mat[i,2])
#
#   idx <- which(aliasSymbol$alias_symbol %in% new_mat[i,2])[1]
#   aliasSymbol$symbol[idx]
# })
# new_mat[,2] <- syn_vec
