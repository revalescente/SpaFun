# plot_smoothed_Lfun Function to plot the smoothed Lcross object

plot_smoothed_Lfun Function to plot the smoothed Lcross object

## Usage

``` r
plot_smoothed_Lfun(
  smoothed_obj,
  clusters_vec,
  mark_i = "Mark 1",
  mark_j = "Mark 2",
  center_plot = TRUE
)
```

## Arguments

- smoothed_obj:

  A `spafun_fd` object from
  [`funSmooth()`](https://revalescente.github.io/SpaFun/reference/funSmooth.md)
  (or a list with elements `r` and `fd`).

- clusters_vec:

  vector of cluster value for each sample, named after the samples

- mark_i:

  name of the mark value for which the Lcross function was calculated

- mark_j:

  name of the mark value for which the Lcross function was calculated

- center_plot:

  to choose if you want the plot centered (L-r ~ r) or not. If the
  curves were already centered in
  [`funSmooth()`](https://revalescente.github.io/SpaFun/reference/funSmooth.md)
  (`centered = TRUE`), they are not centered again.

## Value

A ggplot object
