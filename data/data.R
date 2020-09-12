#' Partitioning of the regions in the brain
#'
#' @name marques
#' @docType data
#' @format A list with two elements
#' \describe{
#'   \item{cell.info}{a data frame with two columns, cell.name and cell.type}
#'   \item{counts}{count matrix with 5069 rows (cells) and 23556 columns (genes)}
#' }
#' @author Kevin Lin \email{kevinl1@wharton.penn.edu}
#' @source All elements of this list were derived from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75330.
#' The counts matrix was derived from the GSE75330_Marques_et_al_mol_counts2.tab file,
#' while the the cell.info data frame was derived from the GSE75330_series_matrix.txt file.
#' @keywords data
NULL


#' Partitioning of the regions in the brain
#'
#' @name zeisel
#' @docType data
#' @format A list with two elements
#' \describe{
#'   \item{cell.info}{a data frame with ten columns, cell_id and level1class}
#'   \item{counts}{count matrix with 3005 rows (cells) and 19972 columns (genes)}
#' }
#' @author Kevin Lin \email{kevinl1@wharton.penn.edu}
#' @source All elements of this list were derived from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361
#' (or alternatively, http://linnarssonlab.org/cortex/). The version in this eSVD package is one that is lightly cleaned
#' in-house.
#' @keywords data
NULL
