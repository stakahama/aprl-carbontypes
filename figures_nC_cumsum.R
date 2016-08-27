
options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(xtable)
library(ggplot2)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R"))
GGTheme()

## -----------------------------------------------------------------------------

DBind[svoc, sk, sj] <- ReadFile("svoc")
DBind[X, Y, Theta, gamma] <- ReadFile("matrices")
molec.moles <- ReadMicromolm3(FilePath("tseries_aer"))
clabels <- ReadFile("clabels")
meas <- ReadFile("meas")$lambdaC
decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

## Y <- Y[svoc,sk]
## X <- X[svoc,sj]
## Theta <- Theta[sk,sj]
## gamma <- gamma[sj]

## -----------------------------------------------------------------------------

colnames(Y) <- clabels[rownames(Theta)]
rownames(Theta) <- clabels[rownames(Theta)]

## -----------------------------------------------------------------------------

## *** calculate nC measured by various tiers of FGs ***

Ystar <- function(Y, Theta, gamma) {
  empty <- !(rowSums(sweep(Theta, 2, gamma, "*"))>0)
  Y[,empty] <- 0
  Y
}

example <- Slice(molec.moles, decisions$hour)
cmpds <- intersect(names(example), rownames(Y))

tiers <- c(head(meas, 1), Map(setdiff, tail(meas, -1), head(meas, -1)))
meas.full <- c(tiers, list("full"=colnames(Theta)))

nC.s <- ldply(meas.full, function(J.s, Y, Theta, gamma, example, cmpds) {
  J.s <- intersect(J.s, colnames(Theta))
  Y.s <- Ystar(Y, Theta[,J.s,drop=FALSE], gamma[J.s])
  prod. <- example[,cmpds] %*% Y.s[cmpds,]
  data.frame(clabel=colnames(Y), nC=c(prod.))
}, Y, Theta, gamma, example, cmpds, .id="tier")

clabel.levs <- with(nC.s %>% filter(tier=="full"),
                    clabel[order(nC, decreasing=TRUE)])

nC.s$clabel <- factor(nC.s$clabel, clabel.levs)

## -----------------------------------------------------------------------------

## *** normalize for plotting ***

wf <- acast(nC.s, clabel~tier, value.var="nC")
wf[] <- apply(wf, 2, cumsum)
wf[,"full"] <- wf[,"full"]-rowSums(wf[,names(meas)])
## colnames(wf) <- c(names(meas), "rest")

wf[] <- wf/sum(wf[nrow(wf),]) # normalize
lf <- melt(wf, c("clabel", "tier"), value.name="cumsum") %>%
  mutate(clabel=factor(clabel, clabel.levs))

levels(lf$tier) <- Capitalize(levels(lf$tier))

ggp <- ggplot(lf %>% filter(clabel %in% clabel.levs[1:20]))+
  geom_bar(aes(clabel, cumsum, fill=tier), stat="identity", position="stack")+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=14))+
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0))+
  labs(x="Carbon type", y="Cumulative carbon fraction")+
  scale_fill_brewer(name="", palette = "Set1")

pdf(FilePath("plot_nC_cumsum"), width=7, height=5.5)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

## *** export table (carbon type abundance and their groups) ***

th <- Theta[clabel.levs[1:15],]
jj <- apply(th > 0, 2, any)
jj <- names(sort(colSums(th[,jj]), decreasing=TRUE))
th.ord <- th[order(as.numeric(rownames(th))),jj]
colnames(th.ord) <- Relabel(colnames(th.ord), labels.FG)
print(xtable(th.ord))

