
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)

File <- function(x, path="data", prefix="merged")
  file.path(path, paste(prefix, x, sep="_"))

## -----------------------------------------------------------------------------

inpfiles <- c(
  "groupattr"=File("group_attributes.csv"),
  "matrices"=File("matrices.rda"),
  "matrices_2"=File("matrices_2.rda"),
  "carbonattr"=File("C_attributes.csv")
)

## -----------------------------------------------------------------------------

load(inpfiles["matrices"])
load(inpfiles["matrices_2"])

carbon.attr <- read.csv(inpfiles["carbonattr"])

## -----------------------------------------------------------------------------

## *** number of functional groups ***

plot(X, sweep(Y %*% Theta, 2, gamma, "*"))
abline(0, 1)

isTRUE(all.equal(X, sweep(Y %*% Theta, 2, gamma, "*")))

## -----------------------------------------------------------------------------

## *** carbon oxidation state per carbon type ***

zeta.table <- unique(subset(carbon.attr,,c("ctype", "OSc")))
zeta.approx <- rowSums(sweep(Theta, 2, gamma*zFG, "*"))

osc <- full_join(zeta.table,
                 data.frame(ctype=names(zeta.approx), zeta=zeta.approx))

with(osc, plot(OSc, zeta))
abline(0, 1)

## which have the discrepancies?
ViewTheta <- function(theta) apply(theta > 0, 1, which)
subset(osc, zeta != OSc)
ViewTheta(Theta[with(osc, ctype[zeta != OSc]),])

## *** mean carbon oxidation state per molecule ***

zeta <- with(zeta.table, setNames(OSc, ctype))[colnames(Y)]

nC <- rowSums(Y)

plot((Y %*% zeta)/nC, (X %*% zFG)/nC)
abline(0,1)

## -----------------------------------------------------------------------------
