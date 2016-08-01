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
  "adjacent"=File("adjacent_atoms.csv")
)

outfiles <- c(
  "groupattr"=File("group_attributes.csv")
)

## -----------------------------------------------------------------------------

fulltable <- read.csv(inpfiles["fulltable"])
adjtable <- read.csv(inpfiles["adjacent"])

fulltable$shorttype <- Shorttype(fulltable$type)
adjtable$atom1_shorttype <- Shorttype(adjtable$atom1_type)
adjtable$atom2_shorttype <- Shorttype(adjtable$atom2_type)

## -----------------------------------------------------------------------------

## atomic composition

numatoms <- fulltable %>%
  count(compound, group, match, shorttype)

numatoms.uniq <- unique(subset(numatoms,,c(group, shorttype, n)))

stopifnot(all(with(numatoms.uniq, table(group, shorttype)) <= 1))

numatoms.wf <- dcast(numatoms.uniq, group ~ shorttype, value.var="n", fill=0)

## -----------------------------------------------------------------------------

## how is it bonded to carbon? based on this, compute oxidation state

ox <- c(
  "C"=0,
  "H"=-1,
  "N"=1,
  "O"=1
)

## match second carbon attached to paired carbon with its
##   associated functional group

bonded <- inner_join(adjtable %>% filter(atom1_shorttype=="C"),
                     fulltable %>% filter(shorttype!="C") %>%
                     rename(atom2=atom) %>% mutate(shorttype=NULL),
                     by=c("compound", "atom2"))

## sum bond contributions

bonded.z <- bonded %>%
  mutate(z=ox[as.character(atom2_shorttype)]*bondorder) %>%
  group_by(compound, atom1, match, group) %>%
  summarize(zFG=sum(z))

## identify unique values for each group

bonded.uniq <- unique(subset(bonded.z,,c("group","zFG")))

## bonded.uniq %>% filter(group=="ester") # there are two values

## this average treatment will be accounted for by gamma

bonded.uniq <- bonded.uniq %>% group_by(group) %>% summarize(zFG=mean(zFG))


## -----------------------------------------------------------------------------

group.attr <- full_join(numatoms.wf, bonded.uniq)

## -----------------------------------------------------------------------------

write.csv(group.attr, outfiles["groupattr"], na="", row.names=FALSE)
