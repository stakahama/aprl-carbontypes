
################################################################################
##
## select_SVOCs.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## see LICENSE_GPLv3.txt
##
################################################################################


options(stringsAsFactors=FALSE)

library(dplyr)
library(Rfunctools)
library(RJSONIO)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", "lib/lib_collapse.R")

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")
molec.attr <- ReadFile("molecattr")

## -----------------------------------------------------------------------------

svoc <- setdiff(with(arrange(molec.attr, logC0), compound[round(logC0) <= 3]),
                Radicalgroups(X))

k <- apply(Y[svoc,] > 0, 2, any)
j <- apply(Theta[k,] > 0, 2, any)

## -----------------------------------------------------------------------------

out <- list(compounds=svoc, ctypes=rownames(Theta)[k], groups=colnames(Theta)[j])

## Y <- Y[svoc,sk]
## X <- X[svoc,sj]
## Theta <- Theta[sk,sj]
## gamma <- gamma[sj]

## -----------------------------------------------------------------------------

cat(toJSON(out), file=FilePath("svoc"))

## -----------------------------------------------------------------------------

## clabels <- with(list(sk=out$ctypes, rest=setdiff(rownames(Theta), out$ctypes)),
##                 setNames(c(sprintf("S%02d",seq_along(sk)),
##                            sprintf("T%02d",seq_along(rest))),
##                          c(sk, rest)))

## cat(toJSON(clabels), file=FilePath("clabels"))
