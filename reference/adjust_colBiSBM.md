# Adjust a colBiSBM on a given point

Adjust a colBiSBM on a given point

## Usage

``` r
adjust_colBiSBM(fitted_bisbmpop, Q, depth = 1L, nb_pass = 1L)
```

## Arguments

- fitted_bisbmpop:

  a fitted bisbmpop, obtained by using the `estimate_colBiSBM`

- Q:

  a vector of size 2, containing the coordinates of the model we want to
  fit

- depth:

  the depth (how far from the center to explore) of the moving window.
  Default to 1.

- nb_pass:

  the number of passes of moving window to perform. Default to 1.

## Value

A bisbmpop object models for the collection of networks. Not the same
object as `fitted_bisbmpop`
