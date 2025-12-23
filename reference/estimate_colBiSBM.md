# Estimate a colBiSBM on a collection of networks

Estimate a colBiSBM on a collection of networks

## Usage

``` r
estimate_colBiSBM(
  netlist,
  colsbm_model,
  net_id = NULL,
  distribution = "bernoulli",
  nb_run = 3L,
  global_opts = list(),
  fit_opts = list(),
  Z_init = NULL,
  sep_BiSBM = NULL
)
```

## Arguments

- netlist:

  A list of matrices.

- colsbm_model:

  Which colBiSBM to use, one of "iid", "pi", "rho", "pirho".

- net_id:

  A vector of string, the name of the networks.

- distribution:

  A string, the emission distribution, either "bernoulli" (the default)
  or "poisson".

- nb_run:

  An integer, the number of run the algorithm do. Default to 3.

- global_opts:

  Global options for the outer algorithm and the output. See details.

- fit_opts:

  Fit options for the VEM algorithm. See details

- Z_init:

  An optional bi-dimensional list of size `Q1_max` x `Q2_max` containing
  for each value a list of two vectors of clusters memberships. Default
  to NULL.

- sep_BiSBM:

  A pre-fitted `sep_BiSBM.` Used to avoid end computations. The best way
  to obtain one is to extract from a fitted bisbmpop object. Defaults to
  NULL.

## Value

A bisbmpop object listing a collection of models for the collection. of
networks

## Details

The list of parameters `global_opts` essentially tunes the exploration
process.

- `nb_cores` integer for number of cores used for parallelization.
  Default is 1

- `verbosity` integer for verbosity (0, 1, 2, 3, 4). Default is 1. 0
  will disable completely the output of the function. Note: you can
  access the \$joint_modelisation_preferred attribute to check which
  modelisation is preferred

- `Q1_max` integer for the max size in row to explore. Default is
  computed with the following formula:
  `floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2)`

- `Q2_max` integer for the max size in columns to explore. Default is
  computed with the following formula:
  `floor(log(sum(sapply(netlist, function(A) ncol(A)))) + 2)`

- `nb_models` the number of models to keep for each values of Q1,Q2.
  Default is 5.

- `depth` specifies how large will the moving window be. Default is 1,
  meaning the window will go from (Q1 - 1, Q2 - 1) to (Q1 + 1, Q2 + 1)
  and all the values in the square defined.

- `plot_details` integer to control the display of the exploration and
  moving window process. Values are 0 or 1. Default is 1.

- `max_pass` the maximum number of moving window passes that will be
  executed. Default is 10.

The list of parameters `fit_opts` are used to tune the Variational
Expectation Maximization algorithm.

- `algo_ve` a string to choose the algorithm to use for the variational
  estimation. Available: "fp"

- `verbosity` an integer to choose the level of verbosity of the fit
  procedure. Defaults to 0. Available: 0,1

- `max_vem_steps` an integer setting the number of Variational
  Expectation-Maximization steps to perform. Defaults to 1000.

- `minibatch` a boolean setting wether to use a "minibatch" like
  approach. If set to TRUE during the VEM the networks will be optimized
  in random orders. If set to FALSE they are optimized in the
  lexicographical order. Default to TRUE.

- `tolerance` a numeric, controlling the tolerance for which a criterion
  is considered converged. Default to 1e-6.

- `greedy_exploration_max_steps` the maximum number of iteration of
  greedy exploration to perform. Defaults to 50.

- `greedy_exploration_max_steps_without_improvement` an integer
  indicating for which number of steps the best model must not change to
  be end greedy exploration. Defaults to 5.

- `kmeans_nstart` an integer indicating the number of random starts to
  use for kmeans in spectral clustering.

- `kmeans_iter_max` an integer indicating the maximum number of
  iterations to use for kmeans in spectral clustering.

## See also

[`clusterize_bipartite_networks()`](https://chabert-liddel.github.io/colSBM/reference/clusterize_bipartite_networks.md),
[`bisbmpop`](https://chabert-liddel.github.io/colSBM/reference/bisbmpop.md),
[`fitBipartiteSBMPop`](https://chabert-liddel.github.io/colSBM/reference/fitBipartiteSBMPop.md),
`browseVignettes("colSBM")`

## Examples

``` r
alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
first_collection <- generate_bipartite_collection(
  nr = 50, nc = 25,
  pi = c(0.5, 0.5), rho = c(0.5, 0.5), alpha = alpha1, M = 2
)
second_collection <- generate_bipartite_collection(
  nr = 50, nc = 25,
  pi = c(0.5, 0.5), rho = c(0.5, 0.5), alpha = alpha2, M = 2
)

netlist <- append(first_collection, second_collection)

if (FALSE) { # \dontrun{
# A collection where joint modelisation makes sense
cl_joint <- estimate_colBiSBM(
  netlist = first_collection,
  colsbm_model = "iid",
  global_opts = list(nb_cores = parallelly::availableCores(omit = 1L))
)
# A collection where joint modelisation doesn't make sense
cl_separated <- estimate_colBiSBM(
  netlist = netlist,
  colsbm_model = "iid",
  global_opts = list(nb_cores = parallelly::availableCores(omit = 1L))
)
} # }
```
