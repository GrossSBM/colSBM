# Perform a spectral clustering

Perform a spectral clustering

## Usage

``` r
spectral_clustering(X, K, kmeans.nstart = 400L, kmeans.iter.max = 50L)
```

## Arguments

- X:

  an adjacency matrix

- K:

  the number of clusters

- kmeans.nstart:

  the number of random starts for the kmeans algorithm. Defaults to 400.
  Ensures consistency of the results.

- kmeans.iter.max:

  the maximum number of iterations for the kmeans algorithm. Defaults to
  50.

## Value

A vector : The clusters labels

A vector : The clusters labels
