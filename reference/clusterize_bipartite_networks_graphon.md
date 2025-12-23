# Clusterize bipartite networks in a bottom up fashion with graphon distance

Clusterize bipartite networks in a bottom up fashion with graphon
distance

## Usage

``` r
clusterize_bipartite_networks_graphon(
  netlist,
  colsbm_model,
  net_id = NULL,
  distribution = "bernoulli",
  nb_run = 3L,
  global_opts = list(),
  fit_opts = list(),
  partition_init = NULL,
  full_inference = FALSE,
  keep_history = FALSE,
  verbose = TRUE,
  temp_save_path = tempfile(fileext = ".Rds")
)
```

## Arguments

- netlist:

  A list of matrices.

- colsbm_model:

  Which colSBM to use, one of "iid", "pi", "rho", "pirho", "delta",
  "deltapi".

- net_id:

  A vector of string, the name of the networks.

- distribution:

  A string, the emission distribution, either "bernoulli" (the default)
  or "poisson"

- nb_run:

  An integer, the number of run the algorithm do. Defaults to 3.

- global_opts:

  Global options for the outer algorithm and the output

- fit_opts:

  Fit options for the VEM algorithm

- partition_init:

  WIP A list of fitted collections from which to start the fusions
  Optional fit init from where initializing the algorithm.

- full_inference:

  The default "FALSE", the algorithm stop once splitting groups of
  networks does not improve the BICL criterion. If "TRUE", then continue
  to split groups until a trivial classification of one network per
  group.

- keep_history:

  A boolean, should the function keep the history of the fusions or not.
  Default to FALSE to reduce file size. Note that the fusion history is
  saved in the temporary file.

- verbose:

  A boolean, should the function be verbose or not. Default to TRUE.

- temp_save_path:

  A string, the path where to save the temporary results. Defaults to a
  temporary file.

## Value

A list of models for the recursive partition of the collection of
networks.
