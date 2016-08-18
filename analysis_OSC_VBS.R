

options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_io.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R", "lib/lib_units.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")

moles.molec <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadTSeries)

clabels <- ReadFile("clabels")

carbon.attr <- ReadFile("carbonattr")

molec.attr <- ReadFile("molecattr")

decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

colnames(Y) <- clabels[rownames(Theta)]
rownames(Theta) <- clabels[rownames(Theta)]

carbon.attr$clabel <- clabels[carbon.attr$ctype]

## -----------------------------------------------------------------------------

examp <- ldply(moles.molec, function(x, i, Y) CarbonProd(Slice(x, i), Y),
               decisions$hour, Y, .id="phase") %>%
  mutate(OC=am["C"]*ppb2micromol(nC), time=NULL)

examp <- left_join(examp, subset(molec.attr,,c(compound, logC0)))
examp <- left_join(examp, subset(carbon.attr,,c(clabel, OSC)))
examp$logC0bin <- round(examp$logC0)

examp <- examp %>% group_by(logC0bin, OSC, phase) %>% summarize(OC=sum(OC)) %>% ungroup()

ctheme <- tail(colorRampPalette(c("darkblue", "lightgray", "darkorange"))(9), -1)

ggplot(examp %>% mutate(OSC=factor(OSC, rev(sort(unique(OSC))))))+
  geom_bar(aes(logC0bin, OC, fill=OSC),
           position="stack", stat="identity")+
  scale_fill_manual(values=ctheme)+
  facet_grid(phase~.)


