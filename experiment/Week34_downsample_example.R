rm(list=ls())
set.seed(10)
load("../../raw_data/marques.RData")

dat <- marques$counts

cell_types <- unique(marques$cell.info$cell.type)
#cell_types <- cell_types[-which(cell_types %in% c("OPC", "PPR"))]

cell_idx <- unlist(lapply(cell_types, function(x){
  tmp <- which(marques$cell.info$cell.type == x)
  sample(tmp, round(length(tmp)/5))
}))
dat <- dat[cell_idx,]
dim(dat)

# use the predetermined set of genes
# gene_mat <- readxl::read_excel("../../raw_data/Marques_genes.xlsx", range = "B3:Y53")
# gene_mat <- as.matrix(gene_mat)
load("../data/marker_genes.Rda")
gene_vec <- as.vector(gene_mat)
gene_vec <- gene_vec[!duplicated(gene_vec)]
gene_idx <- which(colnames(dat) %in% gene_vec)
dat <- dat[,gene_idx]

column_vec <- colnames(dat)[which(colnames(dat) %in% gene_vec)]
gene_vec <- gene_vec[which(gene_vec %in% colnames(dat))]
reordered_idx <- order(column_vec)[rank(gene_vec)]
dat <- dat[,reordered_idx]
dim(dat)

# png("../figure/experiment/Week33_marques.png", height = 1800, width = 1000, res = 300, units = "px")
# par(mar = rep(0.5, 4))
# .plot_singlecell(dat)
# graphics.off()

zz <- apply(dat, 2, function(x){length(which(x!=0))})
dat <- dat[,which(zz > 30)]

#####################

set.seed(10)
dat2 <- singlecell::downsample(dat, downsample_rate = 0.5, dropoff_rate = 0.2)
df2 <- as.data.frame(t(dat2$dat))
res <- VIPER::VIPER(df2, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
                    report = FALSE, outdir = NULL, prefix = NULL)
rownames(df2) <- colnames(dat)


library(VIPER)
data(grun)
dat2 <- singlecell::downsample(as.matrix(t(gene.expression)), downsample_rate = 0.5, dropoff_rate = 0.2)
df2 <- as.data.frame(t(dat2$dat))
res <- VIPER::VIPER(df2, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
                    report = FALSE, outdir = NULL, prefix = NULL)

### save the examples
df <- as.data.frame(t(dat))

write.csv(df, file = "../experiment/singlecell.csv")
write.csv(df2, file = "../experiment/singlecell_downsample.csv")

#####################

df_b <- read.csv("https://raw.githubusercontent.com/linnylin92/VIPER/master/data/singlecell.csv", row.names = 1)
df2_b <- read.csv("https://raw.githubusercontent.com/linnylin92/VIPER/master/data/singlecell_downsample.csv", row.names = 1)
res <- VIPER::VIPER(df2, num = 5000, percentage.cutoff = 0.5, minbool = FALSE, alpha = 1,
                    report = FALSE, outdir = NULL, prefix = NULL)
