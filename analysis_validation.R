
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)

## -----------------------------------------------------------------------------

load(FilePath("matrices"))
load(FilePath("matrices_2"))

carbon.attr <- ReadFile("carbonattr")

## -----------------------------------------------------------------------------

## *** number of functional groups ***

plot(X, sweep(Y %*% Theta, 2, gamma, "*"))
abline(0, 1)

isTRUE(all.equal(X, sweep(Y %*% Theta, 2, gamma, "*")))

## -----------------------------------------------------------------------------

## *** carbon oxidation state per carbon type ***

zeta.table <- unique(subset(carbon.attr,,c("ctype", "OSC")))
zeta.approx <- rowSums(sweep(Theta, 2, gamma*zFG, "*"))

osc <- full_join(zeta.table,
                 data.frame(ctype=names(zeta.approx), zeta=zeta.approx))

with(osc, plot(OSC, zeta))
abline(0, 1)

## which have the discrepancies?
ViewTheta <- function(theta) apply(theta > 0, 1, which)
subset(osc, zeta != OSC)
ViewTheta(Theta[with(osc, ctype[zeta != OSC]),])

## *** mean carbon oxidation state per molecule ***

zeta <- with(zeta.table, setNames(OSC, ctype))[colnames(Y)]

nC <- rowSums(Y)

plot((Y %*% zeta)/nC, (X %*% zFG)/nC)
abline(0,1)

## -----------------------------------------------------------------------------
