# Check Fit Options

This function checks if `fit_opts` is a list.

## Usage

``` r
check_fit_opts(
  fit_opts,
  arg = rlang::caller_arg(fit_opts),
  call = rlang::caller_env()
)
```

## Arguments

- fit_opts:

  A list of fit options.

- arg:

  The argument name for error messages (default is the name of
  `fit_opts`).

- call:

  The calling environment (default is the caller environment).

## Value

Throws an error if the checks fail.
