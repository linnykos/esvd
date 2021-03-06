% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slingshot_curves.R
\name{.initialize_weight_matrix}
\alias{.initialize_weight_matrix}
\title{Initialize the weight matrix}
\usage{
.initialize_weight_matrix(cluster_mat, lineages)
}
\arguments{
\item{cluster_mat}{a 0-1 matrix with n rows and k columns}

\item{lineages}{output of \code{.get_lineages()}}
}
\value{
0-1 matrix of size n by \code{length(lineages)}
}
\description{
Initialize the weight matrix
}
