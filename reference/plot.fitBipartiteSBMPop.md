# Plot matrix summaries of the collection mesoscale structure

Plot matrix summaries of the collection mesoscale structure

## Usage

``` r
# S3 method for class 'fitBipartiteSBMPop'
plot(
  x,
  type = "graphon",
  oRow = NULL,
  oCol = NULL,
  mixture = FALSE,
  net_id = 1L,
  values = FALSE,
  values_min = 0.2,
  ...
)
```

## Arguments

- x:

  a fitBipartiteSBMPop object.

- type:

  The type of the plot. Could be "graphon", "meso" or "block".

- oRow:

  A reordering of the row blocks.

- oCol:

  A reordering of the column blocks.

- mixture:

  Should the block proportions of each network be plotted as well?

- net_id:

  Use to plot only on network in "graphon" view.

- values:

  Wether or not to plot values on the alpha, pi and rho representation.

- values_min:

  The minimum numeric value to plot the value on proportions plots

- ...:

  Further argument to be passed

## Value

A plot, a ggplot2 object.

## Examples

``` r
alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)

first_collection <- generate_bipartite_collection(
  nr = 50, nc = 25,
  pi = c(0.5, 0.5), rho = c(0.5, 0.5),
  alpha = alpha1, M = 2
)

if (FALSE) { # \dontrun{
# A collection where joint modelisation makes sense
cl_joint <- estimate_colBiSBM(
  netlist = first_collection,
  colsbm_model = "iid",
  global_opts = list(nb_cores = parallelly::availableCores(omit = 1L))
)
plot(cl_joint$best_fit)
} # }
```
