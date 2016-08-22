
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

stopifnot(isTRUE(all.equal(X, sweep(Y %*% Theta, 2, gamma, "*"))))

## -----------------------------------------------------------------------------

## *** OM ***

mw <- with(molec.attr, setNames(MW, compound))
cOM <- Ctypemass(Theta, gamma, Lambda)

plot(mw[rownames(Y)], Y %*% cOM[colnames(Y)])
abline(0, 1)

all.equal(mw[rownames(Y)], (Y %*% cOM[colnames(Y)])[,1])

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

stopifnot(isTRUE(all.equal((Y %*% zeta)/nC, (X %*% zFG)/nC)))

## -----------------------------------------------------------------------------

J.s <- c("alkane CH", "alcohol", "ester")
C.s <- rowSums(Theta[,J.s]) > 0
nC.1 <- rowSums(Y[,C.s])
nC.2 <- rowSums(sweep(Y, 2, sign(Theta[,j] %*% gamma[j]), "*"))

plot(nC.1, nC.2)
abline(0, 1)

stopifnot(isTRUE(all.equal(nC.1, nC.2)))
