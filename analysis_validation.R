
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", "lib/lib_metrics.R")

## -----------------------------------------------------------------------------

load(FilePath("matrices"))
load(FilePath("matrices_2"))

carbon.attr <- ReadFile("carbonattr")
molec.attr <- ReadFile("molecattr")

## -----------------------------------------------------------------------------

## *** number of functional groups ***

plot(X, sweep(Y %*% Theta, 2, gamma, "*"))
abline(0, 1)

isTRUE(all.equal(X, sweep(Y %*% Theta, 2, gamma, "*")))

## -----------------------------------------------------------------------------

## *** OM ***

mw <- with(molec.attr, setNames(MW, compound))
cOM <- Ctypemass(Theta, gamma, Lambda)

plot(mw[rownames(Y)], Y %*% cOM[colnames(Y)])
abline(0, 1)




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
                                        # ester-containing groups
## does not agree per carbon type because of the way esters are defined now,
##   but is corrected when summed over a molecule

## -----------------------------------------------------------------------------

## *** mean carbon oxidation state per molecule ***

zeta <- with(zeta.table, setNames(OSC, ctype))[colnames(Y)]

nC <- rowSums(Y)

plot((Y %*% zeta)/nC, (X %*% zFG)/nC)
abline(0,1)

isTRUE(all.equal((Y %*% zeta)/nC, (X %*% zFG)/nC))

## -----------------------------------------------------------------------------
