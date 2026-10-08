# fsmooth Function to convert a vector of data points into a continuous functional object

fsmooth Function to convert a vector of data points into a continuous
functional object

## Usage

``` r
fsmooth(
  dt,
  M = 6,
  genLfun.fd = 4,
  centered = FALSE,
  nbasis = NULL,
  breaks = NULL,
  lambda = NULL,
  rangeval = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- dt:

  object of class data.table with dimensions number of data point per
  number of samples

- M:

  spline order (norder)

- genLfun.fd:

  Lfdobj for fdPar

- centered:

  center by r before smoothing

- nbasis:

  allow user to set number of basis functions

- breaks:

  provide a coarser breaks vector

- lambda:

  smoothing penalty; if NULL a default is chosen

- rangeval:

  override range(r) if desired

- verbose:

  verbose true or false

- ...:

  other arguments for smooth.basis function

## Value

A list with the fd objects and the original matrix and r vector
