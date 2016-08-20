
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_collapse.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")
molec.attr <- ReadFile("molecattr")
measlist <- ReadFile("meas")$lambdaC
svoc <- ReadFile("svoc")$compounds

## -----------------------------------------------------------------------------

## (unnecessary)
## collapserule <- fromJSON(inpfiles["meas"])$collapse
## DBind[measlist.collapsed, matrices.collapsed] <-
##   AggGroups(collapserule, measlist, matrices)

## -----------------------------------------------------------------------------

nC <- rowSums(Y[svoc,])

for(.label in names(measlist)) {

  meas <- measlist[[.label]]

  measC <- rowSums(Theta[,meas]) > 0
  Theta.s <- sweep(Theta[measC, meas], 2, gamma[meas],`*`)
  nC.s <- rowSums(sweep(Y[svoc,measC], 2, sign(rowSums(Theta.s)), `*`))
  Y.s <- Y[svoc,measC]
  X.s <- X[svoc,meas]
  gamma.s <- gamma[meas]

  ## -------------------------------------------------------------------------

  ## *** compound nC comparison ***

  paired <- data.frame(nC, nC.s)

  ggp <- ggplot(paired)+
    geom_point(aes(nC, nC.s))+
    lims(x=c(0, 15), y=c(0, 15))+
    labs(x=expression(n[C]), y=expression(n[C]*"*"))+
    geom_abline(intercept=0, slope=1)

  pdf(SprintF("plot_compound_nC", .label))
  print(ggp)
  dev.off()

}
