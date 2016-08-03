
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("mylib", "lib/lib_collapse.R")

## -----------------------------------------------------------------------------

AddPrefix <- function(x, prefix="merged")
  paste(prefix, x, sep="_")

OutFile <- function(x, ext=NULL, path="outputs", prefix="nC")
  function(suffix=NULL)
    file.path(path, paste(sep=".", paste(sep="_", prefix, x, suffix), ext))

inpfiles <- c(
  "molecattr"=file.path("data", AddPrefix("molec_attributes.csv")),
  "matrices"=file.path("data", AddPrefix("matrices.rda")),
  "meas"=file.path("inputs", "meas_FGs.json"),
  "svoc"=file.path("inputs", "SVOCs.json")
)

outfiles <- list(
  "plot_compound_nC"=OutFile("nC", "pdf")
)

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadRDA(inpfiles["matrices"])

molec.attr <- read.csv(inpfiles["molecattr"])

measlist <- fromJSON(inpfiles["meas"])$lambdaC

svoc <- fromJSON(inpfiles["svoc"])$compounds

## -----------------------------------------------------------------------------

## (unnecessary)
## collapserule <- fromJSON(inpfiles["meas"])$collapse
## DBind[measlist.collapsed, matrices.collapsed] <-
##   AggGroups(collapserule, measlist, matrices)

## -----------------------------------------------------------------------------

nC <- rowSums(matrices$Y[svoc,])

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

  pdf(outfiles[["plot_compound_nC"]](.label))
  print(ggp)
  dev.off()

}
