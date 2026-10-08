# Inverse-density (KDE) downsampling

Samples points with probability proportional to `1 / density`, where the
density is a Gaussian kernel estimate. Points in sparse areas are kept
more often, so the result is closer to a homogeneous pattern.

## Usage

``` r
downsample_KDE(
  pp,
  n_target = NULL,
  prop = NULL,
  sigma = spatstat.explore::bw.scott,
  dimyx = 256,
  ...
)
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

- sigma:

  Kernel bandwidth. Either a number or a function that takes a `ppp` and
  returns a bandwidth. Default
  [`spatstat.explore::bw.scott()`](https://rdrr.io/pkg/spatstat.explore/man/bw.scott.html),
  which is fast on large patterns.

- dimyx:

  Resolution of the pixel grid used for the density estimate (passed to
  [`spatstat.explore::density.ppp()`](https://rdrr.io/pkg/spatstat.explore/man/density.ppp.html)).
  The density is computed on the grid and looked up at the points, which
  keeps the cost almost independent of the number of points.

- ...:

  Further arguments passed to
  [`spatstat.explore::density.ppp()`](https://rdrr.io/pkg/spatstat.explore/man/density.ppp.html).

## Value

A `ppp` with attribute `"downsample"`.

## Details

**Warning:** this changes the second-order structure (K, L, pair
correlation) of the pattern and reduces the clustering you would
measure. Use it for visualisation or for methods that need even
coverage, not before estimating `Lcross`; use
[`downsample_uniform()`](https://revalescente.github.io/SpaFun/reference/downsample_uniform.md)
instead.

## Examples

``` r
library(spatstat.random)
set.seed(1)
X <- rMatClust(20, 0.05, 100)
Y <- downsample_KDE(X, prop = 0.2)
spatstat.geom::npoints(Y)
#> [1] 453
```
