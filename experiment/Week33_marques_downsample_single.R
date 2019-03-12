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

png("../figure/experiment/Week33_marques.png", height = 1800, width = 1000, res = 300, units = "px")
par(mar = rep(0.5, 4))
.plot_singlecell(dat)
graphics.off()

zz <- apply(dat, 2, function(x){length(which(x!=0))})
dat <- dat[,which(zz > 30)]

png("../figure/experiment/Week33_marques.png", height = 1800, width = 1000, res = 300, units = "px")
par(mar = rep(0.5, 4))
.plot_singlecell(dat)
graphics.off()

###############

trials <- 10
downsample_list <- vector("list", trials)
viper_list <- vector("list", trials)
our_list <- vector("list", trials)
i <- 1
set.seed(10*i)
downsample_list[[i]] <- singlecell::downsample(dat, downsample_rate = 0.3, dropoff_rate = 0.2)

png("../figure/experiment/Week33_marques_downsample.png", height = 1800, width = 1000, res = 300, units = "px")
par(mar = rep(0.5, 4))
.plot_singlecell(downsample_list[[i]]$dat)
graphics.off()

quantile(apply(downsample_list[[i]]$dat, 2, function(x){length(which(x!=0))}))
quantile(apply(downsample_list[[i]]$dat, 1, function(x){length(which(x!=0))}))

tmp <- as.data.frame(t(downsample_list[[i]]$dat))
viper_list[[i]] <- VIPER(tmp, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
                    report = FALSE, outdir = NULL, prefix = NULL)
#                    selection = "simplified_expectation")
save.image("Week33_downsample.RData")
