
################################################################################
##
## figures_nC_cumsum.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################


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
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R"))
GGTheme()

## -----------------------------------------------------------------------------

DBind[svoc, sk, sj] <- ReadFile("svoc")
DBind[X, Y, Theta, gamma, zFG, Lambda] <- c(ReadFile("matrices"), ReadFile("matrices_2"))
molec.moles <- ReadMicromolm3(FilePath("tseries_aer"))
clabels <- ReadFile("clabels")
meas <- ReadFile("meas")$lambdaC
decisions <- as.list(ReadFile("example_1"))
carbon.attr <- ReadFile("carbonattr")

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

Psi <- Theta %*% t(Lambda[heteroatoms,])
cOM <- (am["C"] + Psi %*% am[heteroatoms])[,1]

## -----------------------------------------------------------------------------

## *** normalize for plotting ***

CumSumNorm <- function(x, id.vars, ..., levs) {
  ##
  wf <- acast(x, paste(collapse="~", id.vars), ...)
  col <- list(head(colnames(wf), -1), tail(colnames(wf), 1))
  wf[] <- apply(wf, 2, cumsum)
  wf[,col[[2]]] <- wf[,col[[2]]]-rowSums(wf[,col[[1]]])
  wf[] <- wf/sum(wf[nrow(wf),]) # normalize
  ##
  lf <- melt(wf, id.vars, value.name="cumsum")
  lf[[id.vars[1]]] <- factor(lf[[id.vars[1]]], levs)
  lf
}

lf <- CumSumNorm(nC.s, id.vars=c("clabel", "tier"), value.var="nC", levs=clabel.levs)
levels(lf$tier) <- Capitalize(levels(lf$tier))

ggp <- ggplot(lf %>% filter(clabel %in% clabel.levs[1:20]))+
  geom_bar(aes(clabel, cumsum, fill=tier), stat="identity", position="stack")+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=14))+
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0))+
  labs(x="Carbon type", y="Cumulative carbon fraction")+
  ## scale_fill_brewer(name="", palette = "Set1")
  scale_fill_manual(name="", values=colors.set)

## pdf(FilePath("plot_nC_cumsum"), width=7, height=5.5)
## print(ggp)
## dev.off()

## -----------------------------------------------------------------------------

if(0)
  with(nC.s %>% filter(tier=="full") %>%
       mutate(clabel=as.numeric(clabel), OM=nC*cOM[as.character(clabel)]) %>%
       arrange(clabel) %>% mutate(nC=cumsum(nC)/sum(nC), OM=cumsum(OM)/sum(OM)),
       matplot(as.numeric(clabel), cbind(nC, OM)))

lf.OM <- CumSumNorm(nC.s %>% mutate(OM=nC*cOM[as.character(clabel)]),
                    id.vars=c("clabel", "tier"), value.var="OM",
                    levs=clabel.levs)
levels(lf.OM$tier) <- Capitalize(levels(lf.OM$tier))

lf <- full_join(lf %>% mutate(variable="OC"), lf.OM %>% mutate(variable="OM"))
lf$variable <- factor(lf$variable)

lett <- with(lf, data.frame(lett=sprintf("%s)", letters[seq(nlevels(variable))]),
                            variable=levels(variable)))

ggp <- ggplot(lf %>% filter(clabel %in% clabel.levs[1:20]))+
  geom_bar(aes(clabel, cumsum, fill=tier), stat="identity", position="stack")+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=14))+
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0))+
  labs(x="Carbon type", y="Cumulative fraction")+
  scale_fill_manual(name="", values=colors.set)+
  facet_grid(variable~.)+
  theme(panel.margin.y=unit(1, "lines"))+
  geom_text(x=-Inf, y=Inf, hjust=-0.1, vjust=1.2, aes(label=lett), data=lett, size=5)

pdf(FilePath("plot_nC_cumsum"), width=7, height=7)
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
