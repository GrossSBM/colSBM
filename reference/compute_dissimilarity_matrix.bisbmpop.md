# Compute the dissimilarity matrix for a collection of bipartite networks

Compute the dissimilarity matrix for a collection of bipartite networks

## Usage

``` r
compute_dissimilarity_matrix.bisbmpop(collection, weight = "max", norm = "L2")
```

## Arguments

- collection:

  A bmpop or bisbmpop object on which to build the dissimilarity matrix

- weight:

  The weighting to apply to the block proportions. One of "max" or
  "mean", defaults to "max".

- norm:

  The norm to use, either one of "L1" or "L2". Defaults to "L2".
