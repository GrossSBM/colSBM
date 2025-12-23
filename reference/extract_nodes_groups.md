# Extract nodes groups from a fitSimpleSBMPop, fitBipartiteSBMPop, bmpop or bisbmpop object

Extract nodes groups from a fitSimpleSBMPop, fitBipartiteSBMPop, bmpop
or bisbmpop object

## Usage

``` r
extract_nodes_groups(fit, arg = rlang::caller_arg(fit))
```

## Arguments

- fit:

  A fitSimpleSBMPop, fitBipartiteSBMPop, bmpop or bisbmpop object

- arg:

  The name of the argument

## Value

A data.frame with columns network, node_name, cluster and node_type
