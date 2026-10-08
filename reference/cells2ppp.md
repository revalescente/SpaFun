# Build a point pattern from a table of cells

Turns a table of cell centroids into a `spatstat` point pattern. The
window is either the bounding box of the cells or a tissue mask, given
as a PNG file (a thumbnail of the slide) or as an `owin`.

## Usage

``` r
cells2ppp(
  cells,
  mask = NULL,
  mask_scale = NULL,
  magnif = NULL,
  thumb_mag = 1.25,
  threshold = 0.5,
  invert = FALSE,
  x_col = "x",
  y_col = "y",
  marks_col = "mark"
)
```

## Arguments

- cells:

  A `data.frame` with the coordinates and, optionally, the labels of the
  cells, e.g. the output of
  [`read_cells_h5ad()`](https://revalescente.github.io/SpaFun/reference/read_cells_h5ad.md).

- mask:

  Optional tissue mask: path to a PNG file or an `owin`. If `NULL` the
  window is the bounding box of the cells.

- mask_scale:

  Number of coordinate units per mask pixel, e.g. `32` for coordinates
  at 40x and a thumbnail at 1.25x. For an `owin`, it is used to rescale
  it (`NULL` means no rescaling).

- magnif, thumb_mag:

  Shortcut for `mask_scale = magnif / thumb_mag`: magnification of the
  coordinates (e.g. 20 or 40) and of the mask thumbnail. Use either
  these or `mask_scale`.

- threshold:

  Pixels with intensity above `threshold` (0-1 scale, averaged over
  colour channels) are tissue.

- invert:

  If `TRUE`, tissue is the pixels *below* `threshold`.

- x_col, y_col, marks_col:

  Column names in `cells`. Set `marks_col` to `NULL` for an unmarked
  pattern.

## Value

A `ppp` object. Cells outside the window or with missing coordinates are
removed (with a message); their number is stored in the attribute
`"n_removed"`.

## Details

Coordinates and mask must use the same orientation: with a PNG the first
pixel row is the top of the image, as with centroids in image
coordinates (y growing downwards).

## Examples

``` r
set.seed(1)
cells <- data.frame(
  x = runif(500, 0, 320), y = runif(500, 0, 320),
  mark = factor(sample(c("tumor", "stroma"), 500, replace = TRUE))
)

# bounding-box window
pp <- cells2ppp(cells)

# tissue mask: 10 x 10 pixel PNG, each pixel = 32 units
m <- matrix(0, 10, 10); m[2:9, 2:9] <- 1
f <- tempfile(fileext = ".png"); png::writePNG(m, f)
pp <- cells2ppp(cells, mask = f, magnif = 40, thumb_mag = 1.25)
#> 188 of 500 cell(s) (38%) outside the window removed.
attr(pp, "n_removed")
#> missing_coords outside_window 
#>              0            188 
```
