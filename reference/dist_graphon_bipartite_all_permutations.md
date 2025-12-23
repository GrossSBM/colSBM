# Graphon distance for bipartite SBM over all permutations of the blocks

Graphon distance for bipartite SBM over all permutations of the blocks

## Usage

``` r
dist_graphon_bipartite_all_permutations(pis, rhos, alphas)
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
the two structures. This function takes into account the possible
permutation of the blocks and returns the minimum distance.
