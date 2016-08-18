
ppb2micromol <- function(xi, p=1, T=298.15) {
  # ppb (\xi) to micromoles per cubic meter
  micro <- 1e6
  R <- 8.206e-5
  micro * 1e-9 * xi * p/(R*T)
}

ppb2microg <- function(xi, mw) {
  ## ppb (\xi) to micorgrams per cubic meter
  mw[names(xi)] * ppb2micromol(xi)
}
