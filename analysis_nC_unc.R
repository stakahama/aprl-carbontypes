## -----------------------------------------------------------------------------

OMOC <- seq(1, 2.5, .1)
OCr <- seq(0, 1.2, .1)
unc <- .1

plot.new()
plot.window(c(1, 2.5), c(1, 2.5))
axis(1)
axis(2)
box()
matlines(OMOC, 1+outer(OMOC-1, (1+unc*c(-1,1)), "/"), type="l", lty=2)

plot.new()
plot.window(c(0, 1.2), c(0, 1.2))
axis(1)
axis(2)
box()
matlines(OCr, outer(OCr, (1+unc*c(-1,1)), "/"), type="l", lty=2)

## -----------------------------------------------------------------------------
