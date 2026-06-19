# Partition of a collection of unipartite networks based on their common mesoscale structures

Partition of a collection of unipartite networks based on their common
mesoscale structures

## Usage

``` r
clusterize_unipartite_networks(
  netlist,
  colsbm_model,
  directed = FALSE,
  net_id = NULL,
  distribution = "bernoulli",
  nb_run = 3L,
  global_opts = list(),
  fit_opts = list(),
  fit_init = NULL,
  full_inference = FALSE,
  verbose = TRUE,
  temp_save_path = tempfile(fileext = ".Rds")
)
```

## Arguments

- netlist:

  A list of matrices.

- colsbm_model:

  Which colSBM to use, one of "iid", "pi", "delta", "deltapi",

- directed:

  A boolean, should the networks be considered as directed or not.

- net_id:

  A vector of string, the name of the networks.

- distribution:

  A string, the emission distribution, either "bernoulli" (the default)
  or "poisson"

- nb_run:

  An integer, the number of run the algorithm do.

- global_opts:

  Global options for the outer algorithm and the output

- fit_opts:

  Fit options for the VEM algorithm

- fit_init:

  Do not use! Optional fit init from where initializing the algorithm.

- full_inference:

  The default "FALSE", the algorithm stop once splitting groups of
  networks does not improve the BICL criterion. If "TRUE", then continue
  to split groups until a trivial classification of one network per
  group.

- verbose:

  A boolean, should the function be verbose or not. Default to TRUE.

- temp_save_path:

  A string, the path where to save the temporary results. Defaults to a
  temporary file.

## Value

A list with four elements:

- partition:

  A list of models giving the best partition.

- cluster:

  A vector of integers giving the cluster of each network.

- elapsed_time:

  The total time taken by the clustering procedure.

- clustering_history:

  A matrix with M columns and has much rows as there are cuts during
  partitioning.

## Details

This function makes call to `estimate_colSBM`.

## See also

[`estimate_colSBM()`](https://chabert-liddel.github.io/colSBM/reference/estimate_colSBM.md),
[`fitSimpleSBMPop`](https://chabert-liddel.github.io/colSBM/reference/fitSimpleSBMPop.md),
`browseVignettes("colSBM")`

## Examples

``` r
alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
first_collection <- generate_unipartite_collection(
  n = 50,
  pi = c(0.5, 0.5),
  alpha = alpha1, M = 2
)
second_collection <- generate_unipartite_collection(
  n = 50,
  pi = c(0.5, 0.5),
  alpha = alpha2, M = 2
)

netlist <- append(first_collection, second_collection)

if (FALSE) { # \dontrun{
cl_separated <- clusterize_networks(
  netlist = netlist,
  colsbm_model = "iid"
)
} # }
```
