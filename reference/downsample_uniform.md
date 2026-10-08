# Uniform random downsampling

Keeps a random subset of points, each point with the same probability.
Independent thinning leaves the K and L functions (including the
cross-type ones) unchanged, so it is the method to use before `Lcross`.

## Usage

``` r
downsample_uniform(pp, n_target = NULL, prop = NULL, exact = TRUE)
```

## Arguments

- pp:

  A `ppp` object, or a list of `ppp` objects.

- n_target:

  Number of points to keep. Mutually exclusive with `prop`. With a list,
  either one value or one value per pattern.

- prop:

  Fraction of points to keep, in (0, 1\]. Mutually exclusive with
  `n_target`.

- exact:

  If `TRUE` (default) keeps exactly `n_target` points (simple random
  sampling without replacement). If `FALSE` uses independent thinning
  ([`spatstat.random::rthin()`](https://rdrr.io/pkg/spatstat.random/man/rthin.html)):
  the number of kept points is random, with expected value `n_target`.

## Value

A `ppp` with attribute `"downsample"`.

## Examples

``` r
library(spatstat.random)
set.seed(1)
X <- rpoispp(2000)
spatstat.geom::npoints(downsample_uniform(X, n_target = 300))
#> [1] 300
```
