% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCluster.R
\name{summary.cluster.df}
\alias{summary.cluster.df}
\title{summary.getCluster}
\usage{
\method{summary}{cluster.df}(getClusterObj)
}
\arguments{
\item{getClusterObj}{}
}
\value{
results from table
}
\description{
summary.getCluster
}
\examples{
X <- demo$X[c(1,2,3), ]
res <- getUpDownCluster(X)
res_cluster <- getCluster(res)
summary.getCluster(res_cluster)
}
