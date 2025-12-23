# Graphon distance for bipartite SBM

Graphon distance for bipartite SBM

## Usage

``` r
graphon_distance_bipartite(pis, rhos, alphas)
```

## Arguments

- pis:

  A list of two probability vectors (row)

- rhos:

  A list of two probability vectors (columns)

- alphas:

  A list of two connectivity matrices

## Value

The graphon distance between two mesoscale structure.

## Details

The graphon distance is computed as the L2 norm between the graphons of
the two structures. Please note that this does not take into account the
possible permutation of the blocks.
