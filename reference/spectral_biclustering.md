# Perform a spectral bi-clustering, clusters by row and by columns independently

Relies on the spectral_clustering function defined above

## Usage

``` r
spectral_biclustering(A, Q, kmeans.nstart = 400L, kmeans.iter.max = 50L)
```

## Arguments

- A:

  a bipartite adjacency matrix

- Q:

  the two numbers of clusters

- kmeans.nstart:

  the number of random starts for the kmeans algorithm. Defaults to 400.
  Ensures consistency of the results.

- kmeans.iter.max:

  the maximum number of iterations for the kmeans algorithm. Defaults to
  50.

## Value

A list of two vectors : The clusters labels. They are accessed using
\$row_clustering and \$col_clustering

A list of two vectors : The clusters labels. They are accessed using
\$row_clustering and \$col_clustering
