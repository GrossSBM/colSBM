% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmpop_distance.R
\name{graphon_distance_unipartite}
\alias{graphon_distance_unipartite}
\title{Graphon distance for unipartite SBM}
\usage{
graphon_distance_unipartite(pis, alphas)
}
\arguments{
\item{pis}{A list of two probability vectors}

\item{alphas}{A list of two connectivity matrices}
}
\value{
The graphon distance between two mesoscale structure.
}
\description{
Graphon distance for unipartite SBM
}
\details{
The graphon distance is computed as the L2 norm between the
graphons of the two structures. Please note that this does not take into
account the possible permutation of the blocks.
}
\keyword{internal}
