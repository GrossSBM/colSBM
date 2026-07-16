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

## Details

The list of parameters `global_opts` essentially tunes the exploration
process.

- `nb_cores` integer for number of cores used for parallelization.
  Default is 1.

- `sbm_init` boolean specifying wether `{sbm}` package should be used to
  initialize the algorithm.

- `spectral_init` boolean specifying wether a spectral clustering should
  be used to initialize the algorithm.

- `nb_init` the number of different inits to perform to start the
  algorithm fit. Default to 10.

- `verbosity` integer for verbosity (0, 1, 2, 3, 4). Default is 1. 0
  will disable completely the output of the function. modelisation is
  preferred

- `Q_max` integer for the max size to explore. Default is computed with
  the following formula:
  `floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2)`

- `nb_models` the number of models to keep for each values of Q. Default
  is 5.

- `depth` specifies how far around the best model the exploration will
  fit models. Default is 3.

- `plot_details` integer to control the display of the exploration
  process. Values are 0 or 1. Default is 1.

- `max_pass` the maximum number of passes that will be executed. Default
  is 10.

- `backend` the parallelization backend to use. Options available are
  "no_mc", "future", "parallel". Default is "future". Note: we plan to
  unify everything behind future.

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
