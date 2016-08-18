
## evaluation of thresholds for simplification

options(stringsAsFactors=FALSE)

library(Rfunctools)
library(RJSONIO)
library(plyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_io.R", "lib/lib_units.R", "lib/lib_metrics.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")

carbon.attr <- ReadFile("carbonattr")
molec.attr <- ReadFile("molecattr")

moles.molec <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadTSeries)

decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

mw <- with(molec.attr, setNames(MW, compound))
common <- intersect(names(moles.molec$aer), molec.attr$compound)
n <- unclass(Slice(moles.molec$aer, decisions$hour))[1,common]
m <- ppb2microg(n, mw)
Yp <- n*Y[common,]

cmass <- with(CarbonTypeMass(carbon.attr), setNames(OM, ctype))

## -----------------------------------------------------------------------------

## *** based on ppb/ppt ***

thresh <- c(0, .001, .01, .1, 1, 10)
counts <- sapply(thresh, function(v, m, n) length(which(n > v)), m, n)
masses <- sapply(thresh, function(v, m, n) sum(m[n > v]), m, n)

par(mfrow=c(2,1))
plot(thresh[-1], counts[-1], log="x")
plot(thresh[-1], masses[-1], log="x")

cbind(thresh, masses/max(masses))

## *** based on microg/m^3 ***

thresh <- c(0, .001, .01, .1, 1, 10)
counts <- sapply(thresh, function(v, m) length(which(m > v)), m)
masses <- sapply(thresh, function(v, m) sum(m[m > v]), m)

par(mfrow=c(2,1))
plot(thresh[-1], counts[-1], log="x")
plot(thresh[-1], masses[-1], log="x")

cbind(thresh, masses/max(masses))

## *** mass of compound ***

n.molec <- sort(n*mw[common], decreasing=TRUE)

thresh <- seq(0, 200, 10)
masses <- sapply(thresh, function(v, n) sum(n[1:v]), n.molec)

par(mfrow=c(1,1))
plot(thresh, masses)

cbind(thresh, masses/max(masses))


## *** based on mass of C only ***

nC <- sort(rowSums(Yp), decreasing=TRUE)

thresh <- seq(0, 200, 10)
masses <- sapply(thresh, function(v, n) sum(n[1:v]), nC)

par(mfrow=c(1,1))
plot(thresh, masses)

cbind(thresh, masses/max(masses))

## *** based on carbon abundance ***

nC <- sort(colSums(Yp), decreasing=TRUE)

plot(nC)

thresh <- seq(0, 60, 5)
masses <- sapply(thresh, function(v, n) sum(n[1:v]), nC)

par(mfrow=c(1,1))
plot(thresh, masses)

cbind(thresh, round(masses/max(masses), 4))

## *** based on carbon type mass abundance ***

mC <- sort(colSums(sweep(Yp, 2, cmass[colnames(Yp)], "*")), decreasing=TRUE)

plot(mC)

thresh <- seq(0, 60, 5)
masses <- sapply(thresh, function(v, n) sum(n[1:v]), mC)

par(mfrow=c(1,1))
plot(thresh, masses)

cbind(thresh, round(masses/max(masses), 4))
