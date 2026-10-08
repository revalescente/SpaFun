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

  list with elements "r", breaks vector, and "fd", object of the fd
  class

- clusters_vec:

  vector of cluster value for each sample

- mark_i:

  name of the mark value for which the Lcross function was calculated

- mark_j:

  name of the mark value for which the Lcross function was calculated

- center_plot:

  to choose if you want the plot centered (L-r ~ r) or not.

## Value

A ggplot object
