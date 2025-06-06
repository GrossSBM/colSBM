% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic-function.R
\name{plot.fitSimpleSBMPop}
\alias{plot.fitSimpleSBMPop}
\title{Plot matrix summaries of the collection mesoscale structure}
\usage{
\method{plot}{fitSimpleSBMPop}(x, type = "graphon", ord = NULL, mixture = FALSE, net_id = 1L, ...)
}
\arguments{
\item{x}{a fitSimpleSBMPOP object.}

\item{type}{The type of the plot. Could be "graphon", "meso" or "block".}

\item{ord}{A reordering of the blocks.}

\item{mixture}{Should the block proportions of each network be plotted as
well?}

\item{net_id}{Use to plot only on network in "graphon" view.}

\item{...}{Further argument to be passed}
}
\value{
A plot, a ggplot2 object.
}
\description{
Plot matrix summaries of the collection mesoscale structure
}
\examples{
# Trivial example with Gnp networks:
Net <- lapply(
  list(.7, .7, .2, .2),
  function(p) {
    A <- matrix(0, 15, 15)
    A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
    A <- A + t(A)
  }
)
\dontrun{
cl <- estimate_colSBM(Net,
  colsbm_model = "delta",
  directed = FALSE,
  distribution = "bernoulli",
  nb_run = 1
)
plot(cl$best_fit)
}
}
