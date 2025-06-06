% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{spectral_clustering}
\alias{spectral_clustering}
\title{Perform a spectral clustering}
\usage{
spectral_clustering(X, K, kmeans.nstart = 400L, kmeans.iter.max = 50L)
}
\arguments{
\item{X}{an adjacency matrix}

\item{K}{the number of clusters}

\item{kmeans.nstart}{the number of random starts for the kmeans algorithm.
Defaults to 400. Ensures consistency of the results.}

\item{kmeans.iter.max}{the maximum number of iterations for the kmeans
algorithm. Defaults to 50.}
}
\value{
A vector : The clusters labels

A vector : The clusters labels
}
\description{
Perform a spectral clustering
}
\keyword{internal}
