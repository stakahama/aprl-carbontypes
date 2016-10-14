
################################################################################
##
## build_matrices_2.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## see LICENSE_GPLv3.txt
##
################################################################################


options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
PopulateEnv("IO", "config_IO.R")

## -----------------------------------------------------------------------------

FillNA <- function(x, value=0) replace(x, is.na(x), value)

## -----------------------------------------------------------------------------

group.attr <- ReadFile("groupattr")
mat <- ReadFile("matrices")

## -----------------------------------------------------------------------------

zFG <- with(group.attr, setNames(FillNA(zFG), group))

Lambda <- t(as.matrix(SetIndex(mutate(group.attr, zFG=NULL), "group")))

gamma2 <- with(group.attr, setNames(1/C, group))

stopifnot(identical(mat$gamma, gamma2))

## -----------------------------------------------------------------------------

save(zFG, Lambda, file=FilePath("matrices_2"))
