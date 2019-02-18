rm(list=ls())

##########
# read in the cell types

con <- file("../../raw_data/Marques/GSE75330_series_matrix.txt")
txt <- readLines(con)
close(con)

length(txt)
zz <- as.numeric(sapply(txt, nchar))
grep("infer", txt)
# grab the cell names
cell_names <- strsplit(txt[30], split = "\\\t")[[1]]
cell_names <- cell_names[-1]
cell_names <- as.vector(sapply(cell_names, function(x){
  substr(x, 2, nchar(x)-1)
}))

# grab the cell types
# substr(txt[43], 1, 5000)
cell_types <- strsplit(txt[43], split = "\\\t")[[1]]
cell_types <- cell_types[-1]
cell_types <- as.vector(sapply(cell_types, function(x){
  substr(x, 22, nchar(x)-1)
}))
table(cell_types) ## VLMC is labeled as PPR

cell.info <- data.frame(cell.name = cell_names, cell.type = cell_types)
cell.info <- cell.info[order(cell_names),]

##########
# read in the cell count data
dat <- read.table("../../raw_data/Marques/GSE75330_Marques_et_al_mol_counts2.tab",
                  stringsAsFactors = F)

cell_names2 <- as.character(dat[1,-1])
gene_names2 <- as.character(dat[-1,1])
dat <- dat[-1,-1]
dat <- t(dat)
dat <- matrix(as.numeric(dat), nrow = nrow(dat), ncol = ncol(dat))
colnames(dat) <- gene_names2
rownames(dat) <- cell_names2
dat <- dat[order(rownames(dat)),]

marques <- list(cell.info = cell.info, counts = dat)
save(marques, file = "../../raw_data/marques.RData")
