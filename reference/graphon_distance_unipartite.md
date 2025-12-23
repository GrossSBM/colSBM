# Graphon distance for unipartite SBM

Graphon distance for unipartite SBM

## Usage

``` r
graphon_distance_unipartite(pis, alphas)
```

## Arguments

- pis:

  A list of two probability vectors

- alphas:

  A list of two connectivity matrices

## Value

The graphon distance between two mesoscale structure.

## Details

The graphon distance is computed as the L2 norm between the graphons of
the two structures. Please note that this does not take into account the
possible permutation of the blocks.
