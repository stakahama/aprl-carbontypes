
ReadTSeries <- function(filename, ...) {
  ts <- read.csv(filename, colClasses="numeric", check.names=FALSE, ...)
  zoo::zoo(ts[,-1], ts[,1])
}

ReadMicromolm3 <- function(...) {
  ts <- ReadTSeries(...)
  zoo::coredata(ts)[] <- ppb2micromolm3(zoo::coredata(ts))
  ts
}

Slice <- function(x, ...)
  UseMethod("Slice")

Slice.zoo <- function(x, i, tol=.Machine$double.eps) {
  ix <- sapply(i, function(x, x0, tol) abs(x-x0) < tol, zoo::index(x), tol)
  x[ix,]
}

ppb2micromolm3 <- function(xi, p=1, T=298.15) {
  # ppb (\xi) to micromoles per cubic meter
  micro <- 1e6
  R <- 8.206e-5
  micro * 1e-9 * xi * p/(R*T)
}

ppb2microgm3 <- function(xi, mw, ...) {
  ## ppb (\xi) to micorgrams per cubic meter
  if(is.null(names(xi))) {
    warning("missing names(xi)")
    names(mw) <- names(xi) <- seq_along(xi)
  }
  mw[names(xi)] * ppb2micromolm3(xi, ...)
}
