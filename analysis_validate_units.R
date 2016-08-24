
options(stringsAsFactors=FALSE)

library(Rfunctools)
PopulateEnv("lib", "lib/lib_units.R")
PopulateEnv("IO", "config_IO.R")

## -----------------------------------------------------------------------------

molec.attr <- ReadFile("molecattr")
molec.ppb <- ReadFile("tseries_aer")

molec.micromolm3 <- ReadMicromolm3(FilePath("tseries_aer"))

## -----------------------------------------------------------------------------

mw <- with(molec.attr, setNames(MW, compound))

common <- intersect(names(molec.ppb), names(mw))

conc <- as.matrix(molec.ppb[, common])
conc[] <- ppb2micromol(conc)

time <- molec.ppb[,"TIME"]
OM <-  conc %*% mw[common]

OM.2 <- molec.micromolm3[,common] %*% mw[common]
plot(time, OM, type="l")
plot(OM, OM.2)
