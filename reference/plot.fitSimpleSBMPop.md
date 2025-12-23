# Plot matrix summaries of the collection mesoscale structure

Plot matrix summaries of the collection mesoscale structure

## Usage

``` r
# S3 method for class 'fitSimpleSBMPop'
plot(x, type = "graphon", ord = NULL, mixture = FALSE, net_id = 1L, ...)
```

## Arguments

- x:

  a fitSimpleSBMPOP object.

- type:

  The type of the plot. Could be "graphon", "meso" or "block".

- ord:

  A reordering of the blocks.

- mixture:

  Should the block proportions of each network be plotted as well?

- net_id:

  Use to plot only on network in "graphon" view.

- ...:

  Further argument to be passed

## Value

A plot, a ggplot2 object.

## Examples

``` r
# Trivial example with Gnp networks:
Net <- lapply(
  list(.7, .7, .2, .2),
  function(p) {
    A <- matrix(0, 15, 15)
    A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
    A <- A + t(A)
  }
)
if (FALSE) { # \dontrun{
cl <- estimate_colSBM(Net,
  colsbm_model = "delta",
  directed = FALSE,
  distribution = "bernoulli",
  nb_run = 1
)
plot(cl$best_fit)
} # }
```
