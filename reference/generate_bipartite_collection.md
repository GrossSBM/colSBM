# Generate collection of bipartite networks

Generate collection of bipartite networks

## Usage

``` r
generate_bipartite_collection(
  nr,
  nc,
  pi,
  rho,
  alpha,
  M,
  model = "iid",
  distribution = "bernoulli",
  return_memberships = FALSE
)
```

## Arguments

- nr:

  the number of row nodes or a vector specifying the number of row nodes
  for each of the M networks

- nc:

  the number of column nodes or a vector specifying the number of column
  nodes for each of the M networks

- pi:

  a vector of probability to belong to the row clusters

- rho:

  a vector of probability to belong to the columns clusters

- alpha:

  the matrix of connectivity between two clusters

- M:

  the number of networks to generate

- model:

  the colBiSBM model to use. Available: "iid", "pi", "rho", "pirho"

- distribution:

  the emission distribution : "bernoulli" or "poisson"

- return_memberships:

  a boolean which choose whether the function returns a list containing
  the memberships and the incidence matrices or just the incidence
  matrices. Defaults to FALSE, only the matrices are returned.

## Value

A list of M lists, which contains : \$incidence_matrix,
\$row_blockmemberships, \$col_blockmemberships

## Details

the model parameters if set to any other than iid will shuffle the
provided pi and rho
