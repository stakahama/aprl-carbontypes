
################################################################################
##
## analysis_ctype_thresh.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## see LICENSE_GPLv3.txt
##
################################################################################


## evaluation of thresholds for simplification

options(stringsAsFactors=FALSE)

library(Rfunctools)
library(RJSONIO)
library(plyr)
library(dplyr)
library(pryr)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_metrics.R", "lib/lib_collapse.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma, zFG, Lambda] <- c(ReadFile("matrices"), ReadFile("matrices_2"))
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
molec.attr <- ReadFile("molecattr")
molec.moles <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadMicromolm3)
decisions <- as.list(ReadFile("example_1"))
clabels <- ReadFile("clabels")
carbon.attr <- ReadFile("carbonattr")

## -----------------------------------------------------------------------------

n <- Slice(molec.moles$aer, decisions$hour)

## -----------------------------------------------------------------------------

## which are the unmeasured?

cmpds <- intersect(colnames(n), rownames(Y))

radicalgroup <- grepl("radical", colnames(Theta))
hasradicalC <- apply(Theta[,radicalgroup] > 0, 1, any)

Y <- Y[,!hasradicalC]
X <- X[,!radicalgroup]
Theta <- Theta[!hasradicalC,!radicalgroup]

## -----------------------------------------------------------------------------

ctypes <- ldply(measlist, function(j, Theta)
  data.frame(ctype=rownames(Theta)[rowSums(Theta[,intersect(j, colnames(Theta))]) > 0]),
  Theta, .id="meas")

cdict <- with(ctypes, split(ctype, meas))
cdict$full <- rownames(Theta)
odict <- with(unique(select(carbon.attr, OSC, ctype)), setNames(OSC, ctype))

table(odict[cdict$set2])
table(odict[cdict$full])

setdiff(cdict$full, cdict$set2)

## -----------------------------------------------------------------------------

jj <- do.call(seq, as.list(unname(sapply(c("^type", "ctype"), grep, names(carbon.attr))+c(1,0))))
groups <- as.matrix(SetIndex(unique(carbon.attr[,jj]), "ctype"))

unmeas <- groups[with(cdict, setdiff(full, set2)),]

ReverseDict <- function(x)
  setNames(names(x), x)

cgr <- apply(unmeas > 0, 1, partial(`[`, colnames(groups)))
cgr <- ReverseDict(cgr)

conc <- (n[,cmpds] %*% Y[cmpds,])[1,]

head(conc)

sum(conc[rownames(unmeas)])/sum(conc)

sum(conc[cgr[c("quaternary carbon", "C non quaternary")]])/sum(conc)
