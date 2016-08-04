library(zoo)

ReadTSeries <- function(filename) {
  ts <- read.csv(filename, colClasses="numeric", check.names=FALSE)
  zoo(ts[,-1], ts[,1])
}

Slice <- function(x, ...)
  UseMethod("Slice")

Slice.zoo <- function(x, i, tol=.Machine$double.eps)
  x[abs(index(x)-i) < tol,]
