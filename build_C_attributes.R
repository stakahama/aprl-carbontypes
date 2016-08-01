
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)
PopulateEnv("mylib", c("lib/lib_C_attributes.R", "lib/lib_OSc.R"))

File <- function(x, path="data", prefix="merged")
  file.path(path, paste(prefix, x, sep="_"))

## -----------------------------------------------------------------------------

inpfiles <- c(
  "fulltable"=File("MCMGroups_atomfulltable.csv"),
  "adjacent"=File("adjacent_atoms.csv"),
  "OSc"=File("Osc_atomfulltable.csv")
)

outfiles <- c(
  "carbonattr"=File("C_attributes.csv")
)

## -----------------------------------------------------------------------------

fulltable <- read.csv(inpfiles["fulltable"])
adjtable <- read.csv(inpfiles["adjacent"])
osctable <- read.csv(inpfiles["OSc"])

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
  AtomOSc(osctable)
)

## -----------------------------------------------------------------------------

write.csv(carbon.attr, outfiles["carbonattr"], row.names=FALSE)
