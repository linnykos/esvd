% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectory_manipulation.R
\name{construct_pseudotime_trajectory_matrix}
\alias{construct_pseudotime_trajectory_matrix}
\title{Construst the pseudotime trajectory matrix}
\usage{
construct_pseudotime_trajectory_matrix(slingshot_res, cluster_labels)
}
\arguments{
\item{slingshot_res}{the output from the function \code{eSVD::slingshot}}

\item{cluster_labels}{vector of cluster labels, where
the cluster labels are consecutive positive integers from 1 to
\code{max(cluster_labels)}}
}
\value{
a dataframe 6 columns, labeled
\code{cell_idx}, \code{pseudotime}, \code{dist_to_curve} (Euclidean distance to the curve,
in the space used for \code{eSVD::slingshot}), \code{consensus} (\code{NA} if the cell is specific
to a trajectory, or \code{TRUE} or \code{FALSE} otherwise, depending on whether or not the
two trajectories agree on the cell's pseudotime), and
\code{status} (\code{"both"} for cells shared between both trajectories, or \code{"1"}
or \code{"2"} depending on which trajectory the cell belongs), and \code{cluster_labels}
}
\description{
Given the results in \code{eSVD::slingshot} (in \code{slingshot_res}), construct a data frame that summarizes
the information, where each cell is its own row. This function is somewhat complicated
since \code{eSVD::slingshot} can reweight clusters (meaning certain cells are over-sampled),
so \code{eSVD:::.extract_pseudotimes} needs to re-map all the results in \code{slingshot_res}
to unique cells.
}
\details{
This function currently works only when \code{slingshot_res} contains exactly two estimated
trajectories.
}
