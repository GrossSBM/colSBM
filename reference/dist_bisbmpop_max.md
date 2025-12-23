# Compute the dissimilarity between 2 mesoscale structures for bipartite SBM

Compute the dissimilarity between 2 mesoscale structures for bipartite
SBM

## Usage

``` r
dist_bisbmpop_max(pi, rho, alpha, delta = c(1, 1), weight = "max", norm = "L2")
```

## Arguments

- pi:

  A list of two probability vectors (row)

- rho:

  A list of two probability vectors (columns)

- alpha:

  A list of two connectivity matrices

- delta:

  A vector of 2 density parameters (optional)

- weight:

  One of "max" (default) or "mean". See details

- norm:

  "L1"or "L2" norm for the computation

## Value

The dissimilarity between two mesoscale structure.

## Details

If weight is "max" then the weight of each block is computed as
`pmax(pi[[1]], pi[[2]])`. If "mean", then we take the average. "max"
penalize to a greater extent the difference in block proportion between
structure.

## Examples

``` r
pi <- list(c(0.5, 0.5), c(0.1, 0.9))
rho <- list(c(0.1, 0.9), c(0.5, 0.5))
alpha <- list(
  matrix(c(
    0.9, 0.1,
    0.1, 0.05
  ), byrow = TRUE, nrow = 2),
  matrix(c(
    0.9, 0.1,
    0.1, 0.05
  ), byrow = TRUE, nrow = 2)
)
dist_bisbmpop_max(pi, rho, alpha)
#> [1] 0
```
