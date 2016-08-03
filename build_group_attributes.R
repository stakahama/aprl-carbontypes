options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_C_attributes.R", "lib/lib_OSc.R"))

## -----------------------------------------------------------------------------

fulltable <- ReadFile("fulltable")
adjtable <- ReadFile("adjacent")

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

write.csv(group.attr, FilePath("groupattr"), na="", row.names=FALSE)
