# data-raw

Scripts and notes used to build the example datasets shipped with SpaFun.
The raw files are **not** tracked in git and are excluded from the package
build (`^data-raw$` in `.Rbuildignore`, `data-raw/TCGA/` in `.gitignore`).

## Source data

Cell segmentations of five TCGA diagnostic whole-slide images (H&E),
obtained with HoVer-Net. They were produced by Ilaria Billato as part of 
her PhD work:

> TODO: add bibliographic reference.

The original slide images are from The Cancer Genome Atlas (TCGA) and are
available through the NCI Genomic Data Commons (https://portal.gdc.cancer.gov/).

| Slide |
|---|
| TCGA-05-4382-01Z-00-DX1.76b49a4c-dbbb-48b0-b677-6d3037e5ce88 |
| TCGA-CF-A8HY-01Z-00-DX1.53F599FE-3A6D-4BAC-92B2-CBE3F183BB71 |
| TCGA-EM-A1YA-01Z-00-DX1.6CEACBCA-9D05-4BF4-BFF2-725F869ABFA8 |
| TCGA-N5-A4RU-01Z-00-DX1.BEDD8FC8-AA33-42A5-A8BB-3B2E078D74B4 |
| TCGA-OU-A5PI-01Z-00-DX2.F424E74E-2137-493A-8C33-ACAEA2969448 |

Raw files (~285 MB) are kept on the HPC only, under:

```
data-raw/TCGA/
├── h5ad/   # one <slide>.h5ad.gz per slide
└── mask/   # one <slide>.png tissue mask per slide
```

TODO: add the full HPC path of the original copy.

### File format

- **h5ad**: AnnData files. Despite the `.gz` extension they are plain HDF5
  (not compressed); `read_cells_h5ad()` detects compression from the magic
  bytes, not from the extension.
  - `obsm/spatial`: cell centroid coordinates (pixels at full resolution, 40x);
  - `obs/type`: categorical cell type, 6 levels: `benign epithelial`,
    `inflammatory`, `necrotic`, `neoplastic`, `no label`, `stromal`;
  - `uns/spatial/HoVer_Net/bbox`: bounding box of the segmented area.
- **mask**: PNG tissue masks at thumbnail magnification (1.25x), so one
  mask pixel = 40 / 1.25 = 32 full-resolution pixels. PNG rows follow image
  coordinates (y increases downwards); no transposition is needed. With
  this scale and orientation fewer than 0.05% of cells fall outside the
  mask on all five slides.

## Derived datasets

Built by `tcga_examples.R`:

- `data/tcga_crops.rda` (~580 KB): a list of five `ppp` objects, one square
  crop per slide (12–30 thousand cells each, all 6 cell types present).
  Each crop lies inside the tissue mask and is placed where `neoplastic`
  and `stromal` cells are most balanced. Documented in `man/tcga_crops.Rd`.
- `inst/extdata/tcga_tiny.h5ad` and `inst/extdata/tcga_tiny_mask.png`: a
  very small h5ad file and its mask, used in the examples of the reading
  functions (`read_cells_h5ad()`, `cells2ppp()`, `ann2ppp()`).

## Rebuilding

From the package root, with the raw files in `data-raw/TCGA/`:

```r
source("data-raw/tcga_examples.R")
```
