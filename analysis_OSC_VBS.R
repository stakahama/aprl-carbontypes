

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
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R", "lib/lib_units.R"))

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

examp <- ldply(moles.molec, function(x, i, Y) CarbonProd(Slice(x, i), Y),
               decisions$hour, Y, .id="phase") %>%
  mutate(OC=am["C"]*nC, time=NULL)

examp <- left_join(examp, subset(molec.attr,,c(compound, logC0)))
examp <- left_join(examp, subset(carbon.attr,,c(clabel, OSC)))
examp$logC0bin <- round(examp$logC0)

examp <- examp %>% group_by(logC0bin, OSC, phase) %>% summarize(OC=sum(OC)) %>% ungroup()

ctheme <- setNames(tail(colorRampPalette(c("darkorange", "lightgray", "darkblue"))(9), -1),
                   unique(examp$OSC))

ggplot(examp %>% mutate(OSC=factor(OSC, rev(sort(unique(OSC))))))+
  geom_bar(aes(logC0bin, OC, fill=OSC),
           position="stack", stat="identity")+
  scale_fill_manual(values=ctheme)+
  facet_grid(phase~.)

levs <- with(examp, c("gas", rev(sort(unique(OSC)))))

aer <- full_join(examp %>% filter(phase=="aer") %>%
                 mutate(OSC=factor(OSC, levs),
                        phase=NULL),
                 examp %>% filter(phase=="gas") %>%
                 mutate(OSC=factor("gas", levs), phase=NULL) %>%
                 group_by(logC0bin, OSC) %>%
                 summarize(OC=sum(OC)))

ggp <- ggplot(aer)+
  geom_bar(aes(logC0bin, OC, fill=OSC), color="black",
           position="stack", stat="identity")+
  scale_fill_manual(values=c("white", ctheme))

pdf("outputs/OSC_VBS.pdf", width=8, height=5)
print(ggp)
dev.off()

ggplot(aer)+
  geom_bar(aes(logC0bin, OC, fill=OSC), color="black",
           position="stack", stat="identity")+
  scale_fill_manual(values=c("white", ctheme))+
  lims(y=c(0, 30000))


## examp$OSC <- with(examp, factor(OSC, rev(sort(unique(OSC)))))

total <- colSums(acast(examp, OSC~logC0bin, sum, value.var="OC"))
## sc <- 1/sum(total)
sc <- 1
arrays <- dlply(examp, .(phase), function(x) acast(x, OSC~logC0bin, value.var="OC"))

## flipc <- function(x) x[,rev(seq(ncol(x))),drop=FALSE]
## flipr <- function(x) x[rev(seq(nrow(x))),,drop=FALSE]

pdf("outputs/OSC_VBS.pdf", width=8, height=5)
par(mfrow=c(1,1))
barplot(sc*total, col="white")
barplot(sc*arrays$aer, col=ctheme[rownames(arrays$aer)], add=TRUE)
legend("topright", rev(names(ctheme)), fill=rev(ctheme), title=expression(OS[C]))
box()
dev.off()
