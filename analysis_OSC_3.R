
options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(ggplot2)
library(gridExtra)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_io.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R", "lib/lib_units.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")

DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]

DBind[,,merged.osc] <- ReadFile(SprintF("props_file", "actual"))


moles.molec <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadTSeries)

clabels <- ReadFile("clabels")

carbon.attr <- ReadFile("carbonattr")

molec.attr <- ReadFile("molecattr")

decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

examp <- ldply(moles.molec, function(X, i, Y)
  CarbonProd(Slice(X, i), Y) %>% mutate(time=NULL),
  decisions$hour, Y, .id="phase")

examp <- left_join(examp,
                   molec.attr %>% mutate(logC0bin=round(logC0)) %>%
                   select(compound, logC0bin), by="compound")

nC <- examp %>% group_by(phase, logC0bin, ctype) %>% summarize(nC=sum(nC)) %>%
  full_join(., carbon.attr %>% select(ctype, OSC), by="ctype")

ggplot(nC)+
  geom_bar(aes(OSC, nC), stat="identity")+
  facet_grid(logC0bin~phase, scale="free_y")

## -----------------------------------------------------------------------------

ctypes <- ldply(measlist, function(j, Theta) data.frame(ctype=rownames(Theta)[rowSums(Theta[,j]) > 0]), Theta, .id="meas")

wf <- full_join(left_join(ctypes, examp),
                examp %>% mutate(meas="full")) %>%
  dcast(., phase+logC0bin+ctype~meas, sum, value.var="nC")

meas. <- c(names(measlist), "full")

wf[,meas.] <- cbind(wf[,meas.[1],drop=FALSE], t(apply(wf[,meas.], 1, diff)))

## -----------------------------------------------------------------------------

nC <- full_join(melt(wf, measure.vars=meas., variable.name="meas", value.name="nC"),
                carbon.attr %>% select(ctype, OSC), by="ctype")

ggplot(nC)+
  geom_bar(aes(OSC, nC, fill=meas), stat="identity", position="stack")+
  facet_grid(logC0bin~phase, scale="free_y")

## -----------------------------------------------------------------------------

nC.agg <- nC %>% group_by(meas, OSC) %>% summarize(nC=sum(nC))

osc <- merged.osc %>%
  ungroup() %>% filter(case=="ideal") %>% mutate(case=NULL) %>%
  filter((meas=="full" & method=="true") | (meas!="full" & method=="approx")) %>%
  mutate(method=NULL)

ggp1 <- ggplot(nC %>% group_by(OSC) %>% arrange(meas))+
  geom_bar(aes(OSC, nC, fill=meas), stat="identity", position="stack")+
  coord_flip()

ggp2 <- ggplot(osc)+
  geom_hline(yintercept=0)+
  geom_segment(y=0, aes(x=meas, xend=meas, yend=value))+
  geom_point(aes(meas, value))+
  lims(y=c(-4, 3))


pdf("outputs/OSC_distr_meas.pdf", width=10, height=6)
grid.arrange(ggp1, ggp2, ncol=2)
dev.off()
