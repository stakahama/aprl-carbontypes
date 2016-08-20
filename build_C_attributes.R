
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
  rename(atom=atom1) %>% mutate(atom2_type=Shorttype(atom2_type)) %>%
  dcast(., compound+atom~atom2_type, length, value.var="atom2")

merged <- full_join(wf.groups, atoms)

## -----------------------------------------------------------------------------

merged$ctype2 <- apply(merged[,names(am)], 1, Int2Key)

atypes <- unique(merged %>% select(type, ctype, ctype2))
dups <- with(atypes, ctype[duplicated(ctype)])
stopifnot(length(dups)==3) # ester, peroxide, ...

merged$ctype2 <- NULL

## -----------------------------------------------------------------------------

id.vars <- c("compound", "atom", "type", "ctype")

carbon.attr <- full_join(
  merged[, c(id.vars, names(am))],
  AtomOSC(adjtable)
)

## -----------------------------------------------------------------------------

write.csv(carbon.attr, FilePath("carbonattr"), row.names=FALSE)
