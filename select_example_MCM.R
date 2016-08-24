
options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(ggplot2)
library(grid)
## theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R"))
GGTheme()

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")
molec.moles <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadMicromolm3)
clabels <- ReadFile("clabels")
carbon.attr <- ReadFile("carbonattr")
molec.attr <- ReadFile("molecattr")
decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

colnames(Y) <- clabels[rownames(Theta)]
rownames(Theta) <- clabels[rownames(Theta)]

carbon.attr$clabel <- clabels[carbon.attr$ctype]

## -----------------------------------------------------------------------------

moles.lf <- ldply(molec.moles, CarbonInnerProd, Y, .id="phase")

## -----------------------------------------------------------------------------

frac.lf <- moles.lf %>% group_by(phase, time) %>%
  mutate(frac=nC/sum(nC), clabel=factor(clabel, seq(nlevels(clabel)))) %>%
  arrange(as.numeric(clabel))

## cumul <- frac.lf %>% filter(phase=="aer") %>%
##   group_by(clabel) %>%
##   summarize(nC=sum(nC)) %>%
##   arrange(desc(nC))

levs <- c(Gas="gas", Aerosol="aer")
frac.lf$phase <- with(frac.lf, factor(phase, levs, names(levs)))

lett <- data.frame(letter=sprintf("%s)", letters[1:2]),
                   phase=factor(names(levs), names(levs)))

ggp <- ggplot(frac.lf)+
  geom_area(aes(time, frac, fill=clabel), size=.1, color="white")+
  facet_grid(phase~.)+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(expand=c(0, 0))+
  labs(x="Hour", y="Carbon fraction")+
  geom_text(aes(x=Inf, y=Inf, label=letter), hjust=1, vjust=1, size=5, data=lett)+
  theme(panel.margin = unit(.8, "lines"))+
  scale_fill_discrete(name="Carbon type")

pdf(FilePath("plot_ctype_tseries"), width=7, height=5.5)
print(ggp)
dev.off()


## -----------------------------------------------------------------------------

example.aer <- Slice(molec.moles$aer, decisions$hour)
example.aer <- CarbonProd(example.aer, Y) %>% mutate(ctype=factor(ctype))
example.aer <- OrderSlice(example.aer) %>% rename(clabel=ctype)

example.subset <- example.aer %>%
  filter(unclass(compound) <= decisions$ncompounds)

cumul <- example.subset %>% group_by(clabel) %>% summarize(nC=sum(nC)) %>%
  ungroup() %>% arrange(desc(nC)) %>% mutate(cfrac=cumsum(nC)/sum(nC))

clabel.keep <- with(cumul, clabel[1:decisions$ncarbon])

ymax <- ceiling(max(example.subset %>% group_by(compound) %>% summarize(nC=sum(nC)) %>% .[["nC"]]))

ggp <- ggplot(example.subset %>%
              filter(clabel %in% clabel.keep) %>%
              arrange(as.numeric(clabel))) +
  geom_bar(aes(compound, nC, fill=clabel), size=.1, col="white", stat="identity")+
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  labs(x="", y=expression(italic(n)[C]~(mu*mole/m^3)))+
  scale_y_continuous(limits=c(0, ymax+.5), expand=c(0, 0))+
  scale_fill_discrete("Carbon\ntype")

pdf(FilePath("plot_compound_abundance"), width=7, height=5.5)
print(ggp)
dev.off()

stop("*** end of script ***")

## -----------------------------------------------------------------------------

osc.lf <- merge(moles.lf, carbon.attr %>% distinct(clabel, OSC), by="clabel") %>%
  group_by(phase, time, OSC) %>% summarize(nC=sum(nC)) %>% ungroup() %>%
  group_by(phase, time) %>% mutate(frac=nC/sum(nC)) %>% ungroup()

## ggp <- ggplot(osc.lf %>% filter(time %in% seq(0, 25, 1)))+
##   geom_bar(aes(time, frac, fill=OSC), stat="identity")+
##   facet_grid(phase~.)+
##   scale_x_continuous(expand=c(0, 0))+
##   scale_y_continuous(expand=c(0, 0))
## print(ggp)

levs <- sort(unique(osc.lf$OSC), decreasing=TRUE)
colorscale <- colorRampPalette(rev(c("#132B43", "#56B1F7")))(length(levs))

ggp <- ggplot(osc.lf %>% mutate(OSC=factor(OSC, levs)))+
  geom_area(aes(time, frac, fill=OSC), color="white", size=0.1)+
  scale_fill_manual(values=colorscale)+
  facet_grid(phase~.)+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(expand=c(0, 0))

pdf(FilePath("plot_OSC_tseries"), width=10, height=7)
print(ggp)
dev.off()
