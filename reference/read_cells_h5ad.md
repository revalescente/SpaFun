# Read cell centroids and labels from an h5ad file

Reads an AnnData file (`.h5ad`, optionally gzip-compressed `.h5ad.gz`)
and returns one row per cell with its coordinates and its label (for
example the cell type). Requires the Bioconductor package `anndataR`.

## Usage

``` r
read_cells_h5ad(
  data_path,
  marks_col = "type",
  spatial_key = "spatial",
  coords_cols = NULL
)
```

## Arguments

- data_path:

  Path to the `.h5ad` or `.h5ad.gz` file.

- marks_col:

  Name of the `obs` column with the cell labels, or `NULL` to read
  coordinates only.

- spatial_key:

  Name of the `obsm` matrix with the coordinates (first column x, second
  column y). Ignored if `coords_cols` is given.

- coords_cols:

  Optional: names of two `obs` columns holding x and y, to use instead
  of `obsm[[spatial_key]]`.

## Value

A `data.frame` with columns `x`, `y` and, if `marks_col` is not `NULL`,
`mark` (a factor). Row names are the cell ids (`obs_names`).

## Examples

``` r
f <- system.file("extdata", "tcga_tiny.h5ad", package = "SpaFun")
if (nzchar(f) && requireNamespace("anndataR", quietly = TRUE)) {
  cells <- read_cells_h5ad(f, marks_col = "type")
  head(cells)
  table(cells$mark)
}
#> 
#> benign epithelial      inflammatory          necrotic        neoplastic 
#>                33                29                 8              1313 
#>          no label           stromal 
#>                31               113 
```
