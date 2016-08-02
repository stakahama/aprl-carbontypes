library(plyr)
library(dplyr)

Shorttype <- function(x)
  toupper(substring(x, 1, 1))

## AtomOSc <- function(oscfulltable, ...) {
##   ## OSc fulltable

##   ## From Ruggeri and Takahama 2016
##   z <- c(
##     "C-H"=-1,
##     "C-C"=0,
##     "C=C"=0,
##     "C-N"=1,
##     "C-O"=1,
##     "C=O"=2
##   )

##   ## Select carbon atoms and assign Ox numbers
##   fulltable <- fulltable %>%
##     filter(Shorttype(type)=="C") %>%
##     mutate(z=z[as.character(group)])

##   ## Carbon oxidation states
##   fulltable <- fulltable %>%
##     group_by(compound, atom) %>%
##     summarize(OSc=sum(z)) %>%
##     ungroup()

##   ## Return
##   fulltable

## }


AtomOSc <- function(adjtable, ...) {

  ox <- c(
    "C"=0,
    "H"=-1,
    "N"=1,
    "O"=1
  )

  adjtable %>% filter(Shorttype(atom1_type)=="C") %>%
    group_by(compound, atom1) %>%
    summarize(OSc=sum(ox[Shorttype(atom2_type)]*bondorder)) %>%
    ungroup() %>%
    rename(atom=atom1)

}

MolecOSc <- function(table, ...) {

  table <- AtomOSc(table)

  ## Compound mean carbon oxidation states
  table %>%
    group_by(compound) %>%
    summarize(OSc=mean(OSc))

}


GroupOSc <- function(adjtable, fulltable) {

  ox <- c(
    "C"=0,
    "H"=-1,
    "N"=1,
    "O"=1
  )

  adj.C <- adjtable %>%
    mutate(atom1_shorttype=Shorttype(atom1_type),
           atom2_shorttype=Shorttype(atom2_type)) %>%
    filter(atom1_shorttype=="C")

  atoms.nonC <- fulltable %>%
    filter(Shorttype(type)!="C") %>%
    rename(atom2=atom)

  ## match second carbon attached to paired carbon with its
  ##   associated functional group
  grouped <- inner_join(adj.C, atoms.nonC, by=c("compound", "atom2"))

  ## sum bond contributions
  ##   if we did `group_by(compound, atom1, match, group)`,
  ##   `group=="ester"` will have two values
  ## note gamma will be necessary to correct for groups with multiple carbon
  grouped %>%
    group_by(compound, match, group) %>%
    summarize(zFG=sum(ox[as.character(atom2_shorttype)]*bondorder))

}
