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
gene_vec <- sort(unique(as.vector(gene_mat)))
gene_idx <- which(colnames(dat) %in% gene_vec)
dat <- dat[,gene_idx]
dim(dat)

###############

trials <- 10
downsample_list <- vector("list", trials)
viper_list <- vector("list", trials)
our_list <- vector("list", trials)

for(i in 1:trials){
  set.seed(10*i)
  downsample_list[[i]] <- singlecell::downsample(dat)

  tmp <- as.data.frame(t(downsample_list[[i]]$dat))
  viper_list[[i]] <- VIPER::VIPER(tmp, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
                      report = FALSE, outdir = NULL, prefix = NULL)
  save.image("Week33_downsample.RData")

  write.csv(tmp, paste0("Week33_downsample_", i, ".csv"))
  scimpute::scimpute(count_path = paste0("Week33_downsample_", i, ".csv"),
           infile = "csv", outfile = "csv",
           out_dir = paste0("Week33_scimpute_", i),
           drop_thre = 0.5, Kcluster = 6, ncores = 1)

  dat <- log(downsample_list[[i]]$dat+1)
  dropout_mat <- singlecell::dropout(dat)
  zero_mat <- singlecell::find_true_zeros(dropout_mat, num_neighbors = 200)
  idx <- which(is.na(zero_mat))
  our_list[[i]] <- singlecell::scImpute(dat, drop_idx = idx, Kcluster = 6,
                                       verbose = T, weight = 1)
  save.image("Week33_downsample.RData")
}
