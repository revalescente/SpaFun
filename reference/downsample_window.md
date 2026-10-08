# Downsampling by sub-windows

Divides the window into a grid of square tiles and adds randomly chosen
tiles until at least `n_target` points are covered. All points inside
the chosen tiles are kept, and the window of the result is the union of
those tiles intersected with the original window, so edge corrections in
`Lcross` stay correct.

## Usage

``` r
downsample_window(
  pp,
  n_target = NULL,
  prop = NULL,
  tile_size = NULL,
  weighted = FALSE
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

- tile_size:

  Side of the square tiles, in the units of `pp`. Default: the side
  giving about 100 tiles over the bounding box.

- weighted:

  If `TRUE`, tiles are drawn with probability proportional to the number
  of points they contain (fewer empty tiles). Default `FALSE`, i.e.
  uniform over non-empty tiles.

## Value

A `ppp` with attribute `"downsample"`; its element `tiles` is the `tess`
of the chosen tiles.

## Details

Local structure is fully preserved, but the observed area shrinks:
border effects grow and the largest usable distance `r` gets smaller (at
most about `tile_size / 2` for an isolated tile).

## Examples

``` r
library(spatstat.random)
set.seed(1)
X <- rpoispp(2000)
Y <- downsample_window(X, prop = 0.25)
plot(Y)
```
