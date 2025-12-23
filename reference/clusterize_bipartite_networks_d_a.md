# Clusterize bipartite networks with descending and ascending steps

Clusterize bipartite networks with descending and ascending steps

## Usage

``` r
clusterize_bipartite_networks_d_a(
  netlist,
  colsbm_model,
  max_nb_steps = 10L,
  net_id = NULL,
  distribution = "bernoulli",
  nb_run = 3L,
  global_opts = list(),
  fit_opts = list(),
  partition_init = NULL,
  full_collection_init = NULL,
  full_inference = FALSE,
  verbose = TRUE,
  temp_save_path = tempdir()
)
```

## Arguments

- netlist:

  A list of matrices.

- colsbm_model:

  Which colBiSBM to use, one of "iid", "pi", "rho", "pirho",

- max_nb_steps:

  An integer, the maximum number of steps to perform.

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

- partition_init:

  Optional partition list, a list of fitted collections (bisbmpop) from
  which to start partitioning

- full_collection_init:

  Optional full collection, a bisbmpop object containing the fit of all
  the networks

- full_inference:

  The default "FALSE", the algorithm stop once splitting groups of
  networks does not improve the BICL criterion. If "TRUE", then continue
  to split groups until a trivial classification of one network per
  group.

- verbose:

  A boolean, should the function be verbose or not. Default to TRUE.

- temp_save_path:

  A string, the path where to save the temporary results. Defaults to a
  temporary directory.
