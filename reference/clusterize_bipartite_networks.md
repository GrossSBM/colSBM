# Partition of a collection of bipartite networks based on their common mesoscale structures

Partition of a collection of bipartite networks based on their common
mesoscale structures

## Usage

``` r
clusterize_bipartite_networks(
  netlist,
  colsbm_model,
  net_id = NULL,
  distribution = "bernoulli",
  nb_run = 3L,
  global_opts = list(),
  fit_opts = list(),
  partition_init = NULL,
  full_collection_init = NULL,
  full_inference = FALSE,
  method = "single",
  verbose = TRUE,
  temp_save_path = tempfile(fileext = ".Rds")
)
```

## Arguments

- netlist:

  A list of matrices.

- colsbm_model:

  Which colBiSBM to use, one of "iid", "pi", "rho", "pirho",

- net_id:

  A vector of string, the name of the networks.

- distribution:

  A string, the emission distribution, either "bernoulli" (the default)
  or "poisson"

- nb_run:

  An integer, the number of run the algorithm do.

- global_opts:

  Global options for the outer algorithm and the output. See
  `estimate_colBiSBM` for more informations on the elements of the list.

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

- method:

  the agglomeration method to be used. This should be (an unambiguous
  abbreviation of) one of `"ward.D"`, `"ward.D2"`, `"single"`,
  `"complete"`, `"average"` (= ), `"mcquitty"` (= ), `"median"` (= ) or
  `"centroid"` (= ).

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

This functions makes call to `estimate_colBiSBM`.

## See also

[`clusterize_unipartite_networks()`](https://chabert-liddel.github.io/colSBM/reference/clusterize_unipartite_networks.md),
[`estimate_colBiSBM()`](https://chabert-liddel.github.io/colSBM/reference/estimate_colBiSBM.md),
[`fitBipartiteSBMPop`](https://chabert-liddel.github.io/colSBM/reference/fitBipartiteSBMPop.md),
`browseVignettes("colSBM")`

## Examples

``` r
alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
first_collection <- generate_bipartite_collection(
  nr = 50, nc = 25,
  pi = c(0.5, 0.5), rho = c(0.5, 0.5),
  alpha = alpha1, M = 2
)
second_collection <- generate_bipartite_collection(
  nr = 50, nc = 25,
  pi = c(0.5, 0.5), rho = c(0.5, 0.5),
  alpha = alpha2, M = 2
)

netlist <- append(first_collection, second_collection)

if (FALSE) { # \dontrun{
cl_separated <- clusterize_bipartite_networks(
  netlist = netlist,
  colsbm_model = "iid",
  global_opts = list(nb_cores = parallelly::availableCores(omit = 1L))
)
} # }
```
