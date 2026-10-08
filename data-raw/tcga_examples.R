# Example data for SpaFun
#
# Builds, from the 5 TCGA slides in data-raw/TCGA/ (not versioned, ~285 MB):
#   - data/tcga_crops.rda: one square crop per slide, all cells inside it
#     (all 6 cell types), as a named list of marked `ppp`
#   - inst/extdata/: a tiny h5ad (a few thousand cells of one slide) and its
#     tissue mask, used by the examples of read_cells_h5ad() / ann2ppp()
#
# Run from the package root:  source("data-raw/tcga_examples.R")
# Needs: anndataR, spatstat.geom, png, usethis

pkgload::load_all()
library(spatstat.geom)

raw_dir    <- "data-raw/TCGA"
magnif     <- 40       # coordinates at 40x
thumb_mag  <- 1.25     # masks are 1.25x thumbnails -> 32 units per pixel
target_n   <- 20000    # approximate number of cells per crop
min_tissue <- 0.98     # a crop must be (almost) entirely tissue

ids <- sub("\\.h5ad\\.gz$", "", list.files(file.path(raw_dir, "h5ad")))
h5ad_path <- function(id) file.path(raw_dir, "h5ad", paste0(id, ".h5ad.gz"))
mask_path <- function(id) file.path(raw_dir, "mask", paste0(id, ".png"))


# ---- choose a square crop ---------------------------------------------------
# Side chosen so that the crop holds about `target_n` cells, given the mean
# density of the slide. Among the candidate squares on a regular grid that
# lie (almost) entirely in the tissue, keep the one where the two main
# populations (neoplastic and stromal) are most balanced.
choose_crop <- function(pp, target_n, min_tissue = 0.98,
                        balance = c("neoplastic", "stromal")) {
  W <- Window(pp)
  side <- sqrt(target_n / intensity(unmark(pp)))
  bb <- Frame(pp)
  xs <- seq(bb$xrange[1], bb$xrange[2] - side, by = side / 2)
  ys <- seq(bb$yrange[1], bb$yrange[2] - side, by = side / 2)
  cand <- expand.grid(x0 = xs, y0 = ys)

  score <- vapply(seq_len(nrow(cand)), function(k) {
    sq <- owin(cand$x0[k] + c(0, side), cand$y0[k] + c(0, side))
    if (area(intersect.owin(W, sq)) / area(sq) < min_tissue) return(-Inf)
    inside <- inside.owin(pp$x, pp$y, sq)
    n <- table(factor(marks(pp)[inside], levels = levels(marks(pp))))
    min(n[balance])
  }, numeric(1))

  if (all(!is.finite(score))) {
    stop("No candidate square is inside the tissue: lower `target_n` ",
         "or `min_tissue`.")
  }
  best <- which.max(score)
  owin(cand$x0[best] + c(0, side), cand$y0[best] + c(0, side))
}


# ---- 1. crops of the 5 slides -------------------------------------------------
tcga_crops <- lapply(ids, function(id) {
  message("Slide ", id)
  cells <- read_cells_h5ad(h5ad_path(id), marks_col = "type")
  pp <- cells2ppp(cells, mask = mask_path(id), magnif = magnif,
                  thumb_mag = thumb_mag)
  sq <- choose_crop(pp, target_n, min_tissue)
  crop <- pp[sq]
  Window(crop) <- intersect.owin(Window(pp), sq)
  # 2 decimals are more than enough (units are 40x pixels) and compress better
  crop$x <- round(crop$x, 2)
  crop$y <- round(crop$y, 2)
  attr(crop, "n_removed") <- NULL
  crop
})
names(tcga_crops) <- substr(ids, 1, 12)   # TCGA patient barcode
print(sapply(tcga_crops, npoints))

usethis::use_data(tcga_crops, compress = "xz", overwrite = TRUE)


# ---- 2. tiny h5ad + mask in inst/extdata --------------------------------------
# Cells of a small square in the smallest slide, written as a real h5ad so
# that read_cells_h5ad() can be run in examples and tests.
id_small <- ids[which.min(file.size(vapply(ids, h5ad_path, "")))]
ann <- anndataR::read_h5ad(h5ad_path(id_small))
xy  <- as.matrix(ann$obsm$spatial)

pp_small <- tcga_crops[[substr(id_small, 1, 12)]]
c0 <- centroid.owin(Window(pp_small))
half <- 1700                                  # 3000 x 3000 units square (~1500 cells)
keep <- which(abs(xy[, 1] - c0$x) <= half & abs(xy[, 2] - c0$y) <= half)
message(length(keep), " cells in the tiny h5ad")

# obs as read by anndataR can hold columns with S4 classes that write_h5ad()
# cannot write back: rebuild it as a plain data.frame (factors + numerics)
obs_sub <- as.data.frame(ann$obs)[keep, , drop = FALSE]
print(sapply(obs_sub, function(v) class(v)[1]))   # diagnostic
to_plain <- function(v) {
  if (is.factor(v) || is.character(v)) return(factor(as.character(v)))
  out <- tryCatch(as.numeric(v), error = function(e) NULL)
  if (is.null(out)) out <- as.character(v)
  out
}
obs_tiny <- data.frame(lapply(obs_sub, to_plain), check.names = FALSE)
rownames(obs_tiny) <- paste0("cell_", keep)

sp_tiny <- unname(xy[keep, , drop = FALSE])
rownames(sp_tiny) <- rownames(obs_tiny)

# only obs + obsm are needed by read_cells_h5ad(); X is left out
tiny <- anndataR::AnnData(
  obs  = obs_tiny,
  var  = data.frame(row.names = character(0)),
  obsm = list(spatial = sp_tiny)
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
anndataR::write_h5ad(tiny, "inst/extdata/tcga_tiny.h5ad", mode = "w")
file.copy(mask_path(id_small), "inst/extdata/tcga_tiny_mask.png",
          overwrite = TRUE)
file.size(c("inst/extdata/tcga_tiny.h5ad", "inst/extdata/tcga_tiny_mask.png"))
