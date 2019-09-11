rm(list=ls())
load("../experiment/Writeup_interm_geneordering.RData")

# extract the curves
## there are some complications: we had upscaled the dataset, so we need to reverse the process
idx <- our_curves$idx
lambda_long <- our_curves$curves[[1]]$lambda_long
length(idx) == length(lambda_long) # check to see if we have the correct representation

length(which(our_curves$curves[[1]]$lambda_long != 0)) == length(our_curves$curves[[1]]$ord) #... ???
# construct a 3-column matrix:
## lambda_long (which tells which cells are in the curve),
## pos (which we need to fill, NA or the ord values)
## idx (which tells us which cells are duplicated)
mat <- matrix(NA, nrow = length(idx), ncol = 3)
mat[,1] <- our_curves$curves[[1]]$lambda_long
counter <- 1
for(i in 1:nrow(mat)){
  if(mat[i,1] != 0) {
    mat[i,2] <- our_curves$curves[[1]]$ord[counter]
    counter <- counter + 1
  }
}
mat[,3] <- our_curves$idx

# remove cells not in the curve
mat <- mat[apply(mat, 1, function(x){!any(is.na(x))}),]

# remove duplicates
mat <- mat[!duplicated(mat[,3]),]

# append which cell type group each cell is in
mat <- cbind(mat, NA)
for(i in 1:nrow(mat)){
  mat[i,4] <- as.numeric(cell_type_vec[mat[i,3]])
}

# find the common indices we want to keep
intersect_vals <- intersect(our_curves$lineages[[1]], our_curves$lineages[[2]])
mat <- mat[which(mat[,4] %in% intersect_vals),]
dim(mat)

# do a quick sanity check to see if values are roughly in order
sapply(1:length(intersect_vals), function(i){
  idx <- intersect_vals[i]
  summary(mat[which(mat[,4] == idx),2])
}) # something seems wrong? i might want to do this sorting manually, within each cell subtype....


