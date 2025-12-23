# Generate collection of unipartite

Generate collection of unipartite

## Usage

``` r
generate_unipartite_collection(
  n,
  pi,
  alpha,
  M,
  distribution = "bernoulli",
  return_memberships = FALSE
)
```

## Arguments

- n:

  the number of nodes or a vector of the nodes per network

- pi:

  a vector of probability to belong to the clusters

- alpha:

  the matrix of connectivity between two clusters

- M:

  the number of networks to generate

- distribution:

  the emission distribution, either "bernoulli" or "poisson"

- return_memberships:

  Boolean, should return memberships or not.

## Value

A list of M lists, which contains : \$adjacency_matrix, \$clustering

## Details

If n is a single value, this value will be replicated for each of the M
networks. If it is a vector, it must be of size M, specifying the number
of nodes for each network.
