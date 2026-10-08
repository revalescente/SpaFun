# Example crops of five TCGA slides

One square region per slide, with all the cells detected inside it. Each
region was chosen inside the tissue mask, where neoplastic and stromal
cells are most balanced, and holds about 12,000-30,000 cells.

## Usage

``` r
tcga_crops
```

## Format

A named list of 5 marked point patterns (`ppp`), named by TCGA patient
barcode. Coordinates are in pixels of the 40x slide; the marks are a
factor with the cell types `benign epithelial`, `inflammatory`,
`necrotic`, `neoplastic`, `no label` and `stromal`. The window of each
pattern is the square intersected with the tissue mask.

## Source

H&E diagnostic slides from The Cancer Genome Atlas
(<https://portal.gdc.cancer.gov>), cells segmented and classified with
HoVer-Net. Built by `data-raw/tcga_examples.R`.
