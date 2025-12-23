# Check Global Options

This function checks if `global_opts` is a list and if `nb_cores` (if
provided) is an integer greater than a threshold.

## Usage

``` r
check_global_opts(
  global_opts,
  arg = rlang::caller_arg(global_opts),
  call = rlang::caller_env()
)
```

## Arguments

- global_opts:

  A list of global options.

- arg:

  The argument name for error messages (default is the name of
  `global_opts`).

- call:

  The calling environment (default is the caller environment).

## Value

Throws an error if the checks fail.
