
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
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R"))

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
cmpds <- names(example)[example > 1e-3]

tiers <- c(head(meas, 1), Map(setdiff, tail(meas, -1), head(meas, -1)))
meas.full <- c(tiers, list("full"=colnames(Theta)))

nC.s <- ldply(meas.full, function(J.s, Y, Theta, gamma, example, cmpds) {
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
colnames(wf) <- c(names(meas), "rest")

wf[] <- wf/sum(wf[nrow(wf),]) # normalize
lf <- melt(wf, c("clabel", "tier"), value.name="cumsum")

ggp <- ggplot(lf)+
  geom_bar(aes(clabel, cumsum, fill=tier), stat="identity", position="stack")+
  theme(axis.text.x=element_text(angle=60, hjust=1))+
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0))

pdf(FilePath("plot_nC_cumsum"), width=12, height=7)
print(ggp)
dev.off()
