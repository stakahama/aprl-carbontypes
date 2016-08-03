
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(Rfunctools)

source("lib/lib_simpol.R")
simpol <- simpolclass$new()

PopulateEnv("mylib", "lib/lib_OSc.R")

## -----------------------------------------------------------------------------

File <- function(x, path="data", prefix="merged")
  file.path(path, paste(prefix, x, sep="_"))

inpfiles <- c(
  "mass"=file.path("data-raw", "^mcm_.+_mass\\.txt$"),
  "simpol"=File("SIMPOLgroups.csv"),
  ## "OSc"=File("OSc_atomfulltable.csv")
  "adjacent"=File("adjacent_atoms.csv")
)

outfiles <- c(
  "molecattr"=File("molec_attributes.csv")
)

## -----------------------------------------------------------------------------

ReadMassfile <- function(x)
  tryCatch(read.table(x, skip=18, col.names=c("compound", "SMILES", "InChI", "MW")),
           error=function(e) NULL)

## -----------------------------------------------------------------------------

## oscgr <- read.csv(inpfiles["OSc"])
adjacent <- read.csv(inpfiles["adjacent"])
simpgr <- as.matrix(read.csv(inpfiles["simpol"], row.names=1))

massfiles <- list.files(dirname(inpfiles["mass"]), basename(inpfiles["mass"]),
                        full.name=TRUE)
masses <- unique(ldply(massfiles, ReadMassfile))
row.names(masses) <- masses$compound

## -----------------------------------------------------------------------------

temperature <- 298.15
masses[rownames(simpgr),"p0"] <- apply(simpgr, 1, simpol$p0.atm, temperature)

masses <- masses %>%
  mutate(massc=simpol$atm2massc(p0, MW, temperature)) %>%
  mutate(logC0=log10(massc))

## -----------------------------------------------------------------------------

masses <- full_join(masses %>% select(compound, SMILES, MW, logC0),
                   MolecOSc(adjacent))

## -----------------------------------------------------------------------------

write.csv(masses, file=outfiles["molecattr"], row.names=FALSE)
