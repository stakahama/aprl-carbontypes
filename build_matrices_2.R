
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)

File <- function(x, path="data", prefix="merged")
  file.path(path, paste(prefix, x, sep="_"))

## -----------------------------------------------------------------------------

inpfiles <- c(
  "groupattr"=File("group_attributes.csv")
)

outfiles <- c(
  "matrices"=File("matrices_2.rda")
)

## -----------------------------------------------------------------------------

FillNA <- function(x, value=0) replace(x, is.na(x), value)

## -----------------------------------------------------------------------------

group.attr <- read.csv(inpfiles["groupattr"])

## -----------------------------------------------------------------------------

zFG <- with(group.attr, setNames(FillNA(zFG), group))

Lambda <- t(as.matrix(SetIndex(mutate(group.attr, zFG=NULL), "group")))

gamma2 <- with(group.attr, setNames(1/C, group))

## -----------------------------------------------------------------------------

save(zFG, Lambda, gamma2, file=outfiles["matrices"])
