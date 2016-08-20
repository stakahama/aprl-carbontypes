

options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(ggplot2)
## theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_io.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R"))

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

moles.molec[] <- lapply(moles.molec, function(x) {
  coredata(x)[] <- ppb2micromol(coredata(x))
  x
})

moles.lf <- ldply(moles.molec, CarbonInnerProd, Y, .id="phase")

## -----------------------------------------------------------------------------

frac.lf <- moles.lf %>% group_by(phase, time) %>%
  mutate(frac=nC/sum(nC))

## cumul <- frac.lf %>% group_by(phase, time) %>%
##   arrange(desc(frac)) %>% mutate(cumul=cumsum(frac)) %>%
##   slice(with(list(x=which(cumul < 1)), c(x, tail(x,1)+1)))

cumul <- frac.lf %>% filter(phase=="aer") %>%
  group_by(clabel) %>%
  summarize(nC=sum(nC)) %>%
  arrange(desc(nC))

GGTheme()

ggp <- ggplot(frac.lf)+
  geom_area(aes(time, frac, fill=clabel))+
  facet_grid(phase~.)+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(expand=c(0, 0))+
  labs(x="Hour", y="Mass fraction")


theme_update(strip.text=element_text(size=14),
             strip.text.x=element_text(vjust=1),
             strip.text.y=element_text(vjust=.5, angle=90), #vjust=0
             strip.background=element_rect(color=NA, fill=NA, linetype=0),
             axis.text=element_text(size=12),
             axis.text.x=element_text(margin=margin(.5, .5, .5, .5, "lines")),
             axis.text.y=element_text(angle=90, hjust=.5, margin=margin(.5, .5, .5, .5, "lines")),
             axis.ticks.length = unit(-0.3, "lines"),
             panel.border=element_rect(color=1, fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())


pdf(FilePath("plot_ctype_tseries"), width=10, height=7)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

osc.lf <- left_join(moles.lf, carbon.attr %>% select(clabel, OSC), by="clabel") %>%
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
  geom_area(aes(time, frac, fill=OSC), color="white", size=0)+
  scale_fill_manual(values=colorscale)+
  facet_grid(phase~.)+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(expand=c(0, 0))

pdf(FilePath("plot_OSC_tseries"), width=10, height=7)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

example.aer <- Slice(moles.molec$aer, decisions$hour)
example.aer <- CarbonProd(example.aer, Y)
example.aer <- OrderSlice(example.aer)

example.subset <- example.aer %>%
  filter(unclass(compound) < decisions$ncompounds) %>%
  mutate(clabel=factor(clabel))

ggp <- ggplot(example.subset) +
  geom_bar(aes(compound, nC, fill=clabel), stat="identity")+
  theme(axis.text.x=element_text(angle=60, hjust=1))

pdf(FilePath("plot_compound_abundance"), width=10, height=7)
print(ggp)
dev.off()

stop("end script")

## -----------------------------------------------------------------------------
## still exploratory below this point
## -----------------------------------------------------------------------------

## id <- "ctype"#c("type", "ctype")
## atoms <- c("C", "H", "N", "O")
## carbontype.masses <- mutate(unique(carbon.attr[,c(id, atoms)]), C=1)
## carbontype.masses[atoms] <- sweep(carbontype.masses[atoms], 2, am[atoms], "*")
## carbontype.masses$OM <- rowSums(carbontype.masses[atoms])
## carbontype.masses$clabel <- clabels[carbontype.masses$ctype]

carbontype.masses <- CarbonTypeMasses(carbon.attr)

examp <- ldply(moles.molec, function(x, i, Y) CarbonProd(Slice(x, i), Y),
               decisions$hour, Y, .id="phase")
examp$time <- NULL

examp <- left_join(examp, subset(molec.attr,,c(compound, logC0)))
examp <- left_join(examp, subset(carbontype.masses,,c(clabel, OM)))
examp$OM <- with(examp, OM*nC)
examp$logC0bin <- round(examp$logC0)

## examp.clabel <- examp %>% group_by(logC0bin, clabel) %>% summarize(OM=sum(OM))
examp.clabel <- examp %>% group_by(logC0bin, clabel, phase) %>% summarize(OM=sum(OM))
examp.clabel$clabel <- factor(examp.clabel$clabel)
examp.phase <- examp %>% group_by(logC0bin, phase) %>% summarize(OM=sum(OM))

ccol <- rainbow(60)
pcol <- c(aer=rgb(1,1,1,1), gas=rgb(.5,.5,.5,.8))

ggplot()+
  geom_bar(aes(logC0bin, OM, fill=clabel),
           data=examp.clabel,
           position="stack", stat="identity")+
  ## geom_bar(aes(logC0bin, OM, fill=phase), color="black",
  ##          data=examp.phase %>% arrange(-order(phase)),
  ##          position="stack", stat="identity")+
  ## scale_fill_manual(values=c(pcol, ccol))
  ## geom_bar(aes(logC0bin, OM), color="black", fill=NA,
  ##          data=examp.phase %>% arrange(-order(phase)),
  ##          position="stack", stat="identity")+
  facet_grid(phase~.)
