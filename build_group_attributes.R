options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)
PopulateEnv("mylib", c("lib/lib_C_attributes.R", "lib/lib_OSc.R"))

## -----------------------------------------------------------------------------

File <- function(x, path="data", prefix="merged")
  file.path(path, paste(prefix, x, sep="_"))

inpfiles <- c(
  "fulltable"=File("MCMGroups_atomfulltable.csv"),
  "adjacent"=File("adjacent_atoms.csv")
)

outfiles <- c(
  "groupattr"=File("group_attributes.csv")
)

## -----------------------------------------------------------------------------

fulltable <- read.csv(inpfiles["fulltable"])
adjtable <- read.csv(inpfiles["adjacent"])

## -----------------------------------------------------------------------------

## atomic composition

numatoms <- fulltable %>%
  mutate(shorttype=Shorttype(fulltable$type)) %>%
  count(compound, group, match, shorttype)

numatoms.uniq <- unique(subset(numatoms,,c(group, shorttype, n)))

stopifnot(all(with(numatoms.uniq, table(group, shorttype)) <= 1))

numatoms.wf <- dcast(numatoms.uniq, group ~ shorttype, value.var="n", fill=0)

## -----------------------------------------------------------------------------

## how is it bonded to carbon? based on this, compute oxidation state
bonded.z <- GroupOSc(adjtable, fulltable)

## identify unique values for each group
bonded.uniq <- unique(subset(bonded.z,,c("group","zFG")))

## -----------------------------------------------------------------------------

group.attr <- full_join(numatoms.wf, bonded.uniq)

## -----------------------------------------------------------------------------

write.csv(group.attr, outfiles["groupattr"], na="", row.names=FALSE)
