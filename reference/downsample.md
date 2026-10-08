# Downsample a spatial point pattern

Wrapper that reduces the number of points of one or more point patterns
with the chosen method. When `pp` is a list of patterns, samples are
processed in parallel with `future.apply` (set the backend with
[`future::plan()`](https://future.futureverse.org/reference/plan.html));
results are reproducible after
[`set.seed()`](https://rdrr.io/r/base/Random.html).

## Usage

``` r
downsample(
  pp,
  method = c("uniform", "kde", "window"),
  n_target = NULL,
  prop = NULL,
  marks = NULL,
  by_marks = FALSE,
  ...
)
```

## Arguments

- pp:

  A `ppp` object, or a list of `ppp` objects.

- method:

  One of `"uniform"`, `"kde"`, `"window"`.

- n_target:

  Number of points to keep. Mutually exclusive with `prop`. With a list,
  either one value or one value per pattern.

- prop:

  Fraction of points to keep, in (0, 1\]. Mutually exclusive with
  `n_target`.

- marks:

  Optional character vector of cell types (mark levels) to keep; the
  other points are dropped before downsampling.

- by_marks:

  If `TRUE`, each cell type is downsampled separately and `n_target` /
  `prop` apply to each type (e.g. `n_target = 40000` keeps up to 40000
  cells *of each type*). Types with fewer points are kept whole. With
  `method = "uniform"` this is still independent thinning, so cross-type
  K and L functions are unchanged. Not available for
  `method = "window"`.

- ...:

  Further arguments passed to the method function.

## Value

A `ppp` (or a list of `ppp`, with the input names) with attribute
`"downsample"` describing the operation.

## Details

Methods:

- `"uniform"`: random thinning, see
  [`downsample_uniform()`](https://revalescente.github.io/SpaFun/reference/downsample_uniform.md).
  Leaves K and L functions unchanged: the recommended choice before
  computing `Lcross`.

- `"kde"`: inverse-density sampling, see
  [`downsample_KDE()`](https://revalescente.github.io/SpaFun/reference/downsample_KDE.md).
  Makes the pattern more homogeneous: do not use it to estimate K or L.

- `"window"`: keeps all points inside random tiles, see
  [`downsample_window()`](https://revalescente.github.io/SpaFun/reference/downsample_window.md).
  Local structure is intact, observed area shrinks.

## Examples

``` r
library(spatstat.random)
#> Loading required package: spatstat.data
#> Loading required package: spatstat.univar
#> spatstat.univar 3.2-0
#> Loading required package: spatstat.geom
#> spatstat.geom 3.8-3
#> spatstat.random 3.5-2
set.seed(1)
X <- rpoispp(2000)
Y <- downsample(X, method = "uniform", prop = 0.1)
attr(Y, "downsample")
#> $method
#> [1] "uniform"
#> 
#> $n_orig
#> [1] 1971
#> 
#> $n_kept
#> [1] 197
#> 

# several samples, in parallel
# future::plan(future::multisession, workers = 4)
Xs <- list(a = rpoispp(2000), b = rpoispp(3000))
Ys <- downsample(Xs, method = "window", n_target = 500)
sapply(Ys, spatstat.geom::npoints)
#>   a   b 
#> 514 506 
```
