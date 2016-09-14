options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_collapse.R", "lib/lib_constrOptim.R"))

## -----------------------------------------------------------------------------

load(FilePath("matrices"))
carbon.attr <- ReadFile("carbonattr")
molec.attr <- ReadFile("molecattr")
molec.moles <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadMicromolm3)
decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

## examp <- Slice(do.call(`+`, molec.moles), decisions$hour)
examp <- t(sapply(molec.moles, Slice, decisions$hour))["aer",,drop=FALSE]

ix <- intersect(colnames(examp), rownames(Y))
nC <- sort(colSums(examp[,ix] %*% Y[ix,]), decreasing=TRUE)

clabels <- setNames(seq_along(nC), names(nC))
cat(toJSON(clabels), file=FilePath("clabels"))

ix <- intersect(colnames(examp), rownames(X))
fg <- colSums(examp[,ix] %*% X[ix,])

order.FG <- setNames(seq_along(fg), names(fg)[order(fg, decreasing=TRUE)])
cat(toJSON(order.FG), file=FilePath("fgorder"))
