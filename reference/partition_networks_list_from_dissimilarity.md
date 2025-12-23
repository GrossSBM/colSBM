# Partition networks according to the dissimilarity matrix computed with `compute_dissimilarity_matrix`.

Partition networks according to the dissimilarity matrix computed with
`compute_dissimilarity_matrix`.

## Usage

``` r
partition_networks_list_from_dissimilarity(
  networks_list,
  dissimilarity_matrix,
  method = "single",
  nb_groups = 2L
)
```

## Arguments

- networks_list:

  The list of networks to partition

- dissimilarity_matrix:

  The dissimilarity matrix computed with `compute_dissimilarity_matrix`.

- method:

  the agglomeration method to be used. This should be (an unambiguous
  abbreviation of) one of `"ward.D"`, `"ward.D2"`, `"single"`,
  `"complete"`, `"average"` (= ), `"mcquitty"` (= ), `"median"` (= ) or
  `"centroid"` (= ).

- nb_groups:

  An integer, the number of groups. Defaults to 2

## Value

A vector, eventually named according to `networks_list` names.

## Details

This functions partition a provided networks list and outputs back a
vector
