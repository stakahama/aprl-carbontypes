
################################################################################
##
## build_molec_attributes.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## see LICENSE_GPLv3.txt
##
################################################################################


options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(Rfunctools)

source("lib/lib_simpol.R")
simpol <- simpolclass$new()

PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", "lib/lib_OSC.R")

## -----------------------------------------------------------------------------

ReadMassfile <- function(x)
  tryCatch(read.table(x, skip=18, col.names=c("compound", "SMILES", "InChI", "MW")),
           error=function(e) NULL)

## -----------------------------------------------------------------------------

adjacent <- ReadFile("adjacent")
simpgr <- ReadFile("simpol")

massfiles <- list.files(dirname(FilePath("mcmmass")), basename(FilePath("mcmmass")),
                        full.names=TRUE)
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
                   MolecOSC(adjacent))

## -----------------------------------------------------------------------------

write.csv(masses, file=FilePath("molecattr"), row.names=FALSE)
