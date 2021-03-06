% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slingshot_curves.R
\name{.construct_cluster_matrix}
\alias{.construct_cluster_matrix}
\title{Construst cluster matrix from cluster labels}
\usage{
.construct_cluster_matrix(cluster_labels)
}
\arguments{
\item{cluster_labels}{vector of cluster labels, where
the cluster labels are consecutive positive integers from 1 to
\code{max(cluster_labels, na.rm = T)}. Can include \code{NA}}
}
\value{
A 0-1 matrix with \code{length(cluster_labels)} rows
and \code{max(cluster_labels)} columns
}
\description{
Construst cluster matrix from cluster labels
}
