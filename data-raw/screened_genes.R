screened_genes <- read.csv("../data-raw/screened_genes.csv")
colnames(screened_genes) <- "genes"
usethis::use_data(screened_genes, overwrite = T)
