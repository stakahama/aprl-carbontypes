options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_io.R", "lib/lib_metrics.R"))

## -----------------------------------------------------------------------------

matrices <- c(ReadFile("matrices"), ReadFile("matrices_2"))

molec.attr <- ReadFile("molecattr")

DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]

svoc <- ReadFile("svoc")$compounds

clabels <- ReadFile("clabels")

moles.molec <- ReadTSeries(FilePath("tseries_aer"))

decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma, zFG, Lambda] <- matrices

cmpds <- intersect(names(moles.molec), rownames(Y))
n.example <- coredata(Slice(moles.molec[,cmpds], decisions$hour))[1,]

full <- Calculate(n.example, X, Y, Theta, gamma, zFG, Lambda, cmpds)

out <- list()

for(.label in names(measlist)) {

  meas <- measlist[[.label]]

  ctypes <- rownames(Theta)[rowSums(Theta[,meas]) > 0]
  star <- Calculate(n.example, X, Y, Theta, gamma, zFG, Lambda, cmpds, ctypes, meas)

  star$C.recovery <- star$nC.total/full$nC.total
  star$OM.recovery <- star$OM/full$OM

  out[[.label]] <- star

}

## -----------------------------------------------------------------------------
## Exploratory beyond this point
## -----------------------------------------------------------------------------

Extract <- function(x)
  with(x$mix, data.frame(nC=nC.total, OSC, OM, OM.OC, as.list(OM.OC.comp)))

wf <- ldply(out, Extract, .id="meas")
lf <- full_join(melt(Extract(full) %>% mutate(meas="full"), "meas"), melt(wf, "meas"))

atomlist <- c("C", "H", "O", "N")

ggplot(lf %>% filter(!variable %in% atomlist))+
  geom_bar(aes(meas, value), stat="identity")+
  facet_grid(variable~., scale="free_y")

ggplot(lf %>% filter(variable %in% atomlist))+
  geom_bar(aes(meas, value, fill=variable), stat="identity", position="stack")

ggplot(lf %>% filter(!variable %in% "C"))+
  geom_bar(aes(meas, value), stat="identity")+
  facet_grid(variable~., scale="free_y")

## -----------------------------------------------------------------------------
