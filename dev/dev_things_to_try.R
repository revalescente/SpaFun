# curve giocattolo

set.seed(123)

n_t     <- 100                       # punti di tempo (righe)
n_curve <- 20                        # curve (colonne)
tt      <- seq(0, 1, length.out = n_t)

amp   <- runif(n_curve, 0.5, 2)
phase <- runif(n_curve, 0, pi)

# curve vere (senza rumore)
Y_true <- sapply(seq_len(n_curve), function(j)
  amp[j] * sin(2 * pi * tt + phase[j]) + 0.5 * cos(6 * pi * tt))

# dati osservati = curve vere + rumore
df <- Y_true + matrix(rnorm(n_t * n_curve, sd = 0.2), n_t, n_curve)
colnames(df) <- paste0("curva_", seq_len(n_curve))

matplot(tt, df, type = "l", lty = 1, col = adjustcolor(1:n_curve, 0.6),
        xlab = "t", ylab = "y", main = "Dati giocattolo")

#### FDA GCV -------------------------------------------------------------------
library(fda)

lambda_ramsay <- function(Y, tt, norder = 4, Lfd = 2,
                          df_range = c(2.5, 0.5 * length(tt)),
                          n_grid = 50, plot = TRUE) {
  Y <- as.matrix(Y)
  stopifnot(nrow(Y) == length(tt), !is.unsorted(tt, strictly = TRUE))
  
  # 1. basi fisse: un nodo per ogni punto osservato (nbasis = n + norder - 2)
  basis  <- create.bspline.basis(range(tt), norder = norder, breaks = tt)
  Lfdobj <- int2Lfd(Lfd)
  
  fit <- function(loglam, y = Y) smooth.basis(tt, y, fdPar(basis, Lfdobj, 10^loglam))
  
  # i df dipendono solo da lambda e dalle basi, non dai dati: basta una curva
  df_of <- function(loglam) fit(loglam, Y[, 1])$df
  
  # 2. range di lambda ricavato dai gradi di liberta' (df decresce con lambda)
  find_loglam <- function(target)
    uniroot(function(ll) df_of(ll) - target, c(-12, 8), extendInt = "downX")$root
  
  ll_min <- find_loglam(df_range[2])   # lambda piccolo -> df alti (poco liscio)
  ll_max <- find_loglam(df_range[1])   # lambda grande  -> df bassi (quasi retta)
  
  # 3. GCV su griglia nel range trovato
  grid <- seq(ll_min, ll_max, length.out = n_grid)
  tab  <- as.data.frame(t(sapply(grid, function(ll) {
    sm <- fit(ll)
    c(loglam = ll, df = sm$df, gcv = sum(sm$gcv))
  })))
  
  i <- which.min(tab$gcv)
  if (i %in% c(1, n_grid)) warning("Minimo GCV sul bordo del range: allarga df_range")
  
  # 4. raffinamento attorno al minimo
  opt  <- optimize(function(ll) sum(fit(ll)$gcv),
                   c(tab$loglam[max(i - 1, 1)], tab$loglam[min(i + 1, n_grid)]))
  best <- fit(opt$minimum)
  
  if (plot) {
    op <- par(mfrow = c(1, 2)); on.exit(par(op))
    plot(tab$df, tab$gcv, type = "b", xlab = "gradi di liberta'", ylab = "GCV")
    abline(v = best$df, col = 2, lty = 2)
    plot(tab$loglam, tab$gcv, type = "b", xlab = "log10(lambda)", ylab = "GCV")
    abline(v = opt$minimum, col = 2, lty = 2)
  }
  
  list(lambda = 10^opt$minimum, df = best$df, gcv = sum(best$gcv),
       fd = best$fd, basis = basis, table = tab)
}

res <- lambda_ramsay(df, tt)
res$lambda   # lambda ottimo
res$df       # gradi di liberta' equivalenti

Y_smooth <- eval.fd(tt, res$fd)
rmse <- function(a, b) sqrt(mean((a - b)^2))
rmse(Y_smooth, Y_true)
rmse(eval.fd(tt, fit_fda$fd), Y_true)
