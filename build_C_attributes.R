
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_C_attributes.R", "lib/lib_OSC.R"))

## -----------------------------------------------------------------------------

fulltable <- ReadFile("fulltable")
adjtable <- ReadFile("adjacent")

## -----------------------------------------------------------------------------

idvars <- c("compound", "atom")#, "type")

## attached groups
wf.groups <- AddCtypeWide(fulltable)

## attached atoms
wf.atoms <- adjtable %>% filter(Shorttype(atom1_type)=="C") %>%
  rename(atom=atom1) %>%
  mutate(shorttype2=Shorttype(atom2_type)) %>%
  do(dcast(., compound+atom~shorttype2, length, value.var="atom2"))
avars <- sort(setdiff(names(wf.atoms), idvars))
wf.atoms$ctype2 <- apply(wf.atoms[, avars], 1, Int2Key)

## -----------------------------------------------------------------------------

merged <- full_join(wf.groups %>% select_(.dots=c(idvars, "type", "ctype")),
                    wf.atoms %>% select_(.dots=c(idvars, "ctype2")))
atypes <- unique(merged %>% select(type, ctype, ctype2))
dups <- with(atypes, ctype[duplicated(ctype)])
stopifnot(length(dups)==3)

## -----------------------------------------------------------------------------

carbon.attr <- full_join(
  wf.groups %>% select_(.dots=c(idvars, "type", "ctype")),
  wf.atoms %>% select_(.dots=c(setdiff(names(wf.atoms), "ctype2")))
)

carbon.attr <- full_join(
  carbon.attr,
  ## AtomOSC(osctable)
  AtomOSC(adjtable)
)

## -----------------------------------------------------------------------------

write.csv(carbon.attr, FilePath("carbonattr"), row.names=FALSE)
