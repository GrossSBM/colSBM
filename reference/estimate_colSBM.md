# Estimate a colSBM on a collection of networks

Estimate a colSBM on a collection of networks

## Usage

``` r
estimate_colSBM(
  netlist,
  colsbm_model,
  net_id = NULL,
  directed = NULL,
  distribution = "bernoulli",
  fit_sbm = NULL,
  nb_run = 3L,
  global_opts = list(),
  fit_opts = list(),
  Z_init = NULL,
  fit_init = NULL
)
```

## Arguments

- netlist:

  A list of matrices.

- colsbm_model:

  Which colSBM to use, one of "iid", "pi", "delta", "deltapi".

- net_id:

  A vector of string, the name of the networks.

- directed:

  A boolean, are the networks directed or not.

- distribution:

  A string, the emission distribution, either "bernoulli" (the default)
  or "poisson"

- fit_sbm:

  A list of fitted model using the `sbm` package. Use to speed up the
  initialization.

- nb_run:

  An integer, the number of run the algorithm do.

- global_opts:

  Global options for the outer algorithm and the output

- fit_opts:

  Fit options for the VEM algorithm

- Z_init:

  An optional list of cluster memberships for q in 1:Q_max. Default to
  NULL.

- fit_init:

  Do not use! Optional fit init from where initializing the algorithm.

## Value

A bmpop object listing a collection of fitted models for the collection
of networks

## See also

[`clusterize_unipartite_networks()`](https://chabert-liddel.github.io/colSBM/reference/clusterize_unipartite_networks.md),
[`bmpop`](https://chabert-liddel.github.io/colSBM/reference/bmpop.md),
[`fitSimpleSBMPop`](https://chabert-liddel.github.io/colSBM/reference/fitSimpleSBMPop.md),
`browseVignettes("colSBM")`

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
} # }
```
