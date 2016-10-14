
################################################################################
##
## build_C_attributes.R
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
library(pryr)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_C_attributes.R", "lib/lib_OSC.R", "lib/lib_metrics.R"))

## -----------------------------------------------------------------------------

fulltable <- ReadFile("fulltable")
adjtable <- ReadFile("adjacent")
load(FilePath("matrices"))
load(FilePath("matrices_2"))

## -----------------------------------------------------------------------------

wf.groups <- AddCtypeWide(fulltable)

## -----------------------------------------------------------------------------

## attached carbon + heteroatoms
##   don't use this for estimation of mass since some atoms
##   (including heteroatoms) are double-counted

atoms <- adjtable %>% filter(Shorttype(atom1_type)=="C") %>%
  rename(atom=atom1, type=atom1_type) %>%
  mutate(atom2_type=factor(Shorttype(atom2_type), names(am))) %>%
  dcast(., compound+atom+type~atom2_type, length, value.var="atom2")

missing <- setdiff(names(am), names(atoms))
atoms[,missing] <- 0

merged <- full_join(wf.groups, atoms)

## -----------------------------------------------------------------------------

merged$ctype2 <- apply(merged[,names(am)], 1, Int2Key)

atypes <- unique(merged %>% select(type, ctype, ctype2))
dups <- with(atypes, ctype[duplicated(ctype)])
## stopifnot(length(dups)==3) # ester, peroxide, ... when using "merged" table including TMB

merged$ctype2 <- NULL

## -----------------------------------------------------------------------------

id.vars <- c("compound", "atom", "type", "ctype")

carbon.attr <- full_join(
  merged,#[, c(id.vars, names(am))],
  AtomOSC(adjtable)
)

## -----------------------------------------------------------------------------

write.csv(carbon.attr, FilePath("carbonattr"), row.names=FALSE)
