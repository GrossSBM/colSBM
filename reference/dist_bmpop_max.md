# Compute the dissimilarity between 2 mesoscale structures

Compute the dissimilarity between 2 mesoscale structures

## Usage

``` r
dist_bmpop_max(
  pi,
  alpha,
  delta = c(1, 1),
  weight = "max",
  norm = "L2",
  directed
)
```

## Arguments

- pi:

  A list of two probability vectors

- alpha:

  A list of two connectivity matrices

- delta:

  A vector of 2 density parameters (optional)

- weight:

  One of "max" (default) or "mean". See details

- norm:

  "L1"or "L2" norm for the computation

- directed:

  Are the structure of the networks directed?

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
dist_bmpop_max(pi, alpha)
#> [1] 0
```
