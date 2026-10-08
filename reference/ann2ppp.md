# Read an h5ad file as a point pattern

Convenience wrapper:
[`read_cells_h5ad()`](https://revalescente.github.io/SpaFun/reference/read_cells_h5ad.md)
followed by
[`cells2ppp()`](https://revalescente.github.io/SpaFun/reference/cells2ppp.md).

## Usage

``` r
ann2ppp(
  data_path,
  mask = NULL,
  mask_scale = NULL,
  magnif = NULL,
  thumb_mag = 1.25,
  marks_col = "type",
  spatial_key = "spatial",
  coords_cols = NULL,
  ...
)
```

## Arguments

- data_path:

  Path to the `.h5ad` or `.h5ad.gz` file.

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

- marks_col:

  Name of the `obs` column with the cell labels, or `NULL` to read
  coordinates only.

- spatial_key:

  Name of the `obsm` matrix with the coordinates (first column x, second
  column y). Ignored if `coords_cols` is given.

- coords_cols:

  Optional: names of two `obs` columns holding x and y, to use instead
  of `obsm[[spatial_key]]`.

- ...:

  Further arguments passed to
  [`cells2ppp()`](https://revalescente.github.io/SpaFun/reference/cells2ppp.md)
  (`threshold`, `invert`).

## Value

A `ppp` object, see
[`cells2ppp()`](https://revalescente.github.io/SpaFun/reference/cells2ppp.md).

## Examples

``` r
f    <- system.file("extdata", "tcga_tiny.h5ad", package = "SpaFun")
mask <- system.file("extdata", "tcga_tiny_mask.png", package = "SpaFun")
if (nzchar(f) && requireNamespace("anndataR", quietly = TRUE)) {
  # the mask covers the whole slide, the tiny file only a small square of it
  pp <- ann2ppp(f, mask = mask, magnif = 40, thumb_mag = 1.25,
                marks_col = "type")
  pp
}
#> Marked planar point pattern: 1527 points
#> Multitype, with levels = 
#>    benign epithelial inflammatory necrotic neoplastic no label stromal
#> window: binary image mask
#> 2736 x 2896 pixel array (ny, nx)
#> enclosing rectangle: [0, 92672] x [0, 87552] units
```
