
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
## GGTheme()

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma, zFG, Lambda] <- c(ReadFile("matrices"), ReadFile("matrices_2"))
moles.obs <- read.csv("inputs/Sax2005_Table2_APIN.csv", check.names=FALSE, row.names="sample")
moles.sim <- ReadMicromolm3(FilePath("tseries_aer"))
DBind[, measlist, collapserule] <- ReadFile("meas")

## -----------------------------------------------------------------------------

Lambda <- Lambda[c("H", "O", "N"),]
colnames(Lambda) <- Relabel(colnames(Lambda), c(labels.FG, setNames("CO", "ketone"))) # We will use the column of ketone for CO as it dominates the condensed phase over aldehyde (the only difference is in the no. of H). For formal aggregation, use AggGroups() in lib/lib_collapse.R, which assumes equal molar portions the mean.
colnames(Theta) <- colnames(Lambda)
meas.set1 <- Relabel(measlist$set1, c(labels.FG, setNames("CO", "ketone")))

lambdaC <- list(
  old=c("aCH"=.66, "COOH"=1, "CO"=1, "aCOH"=0.66, "CONO2"=0.66),
  new=c("aCH"=.45, "COOH"=1, "CO"=1, "aCOH"=0.5, "CONO2"=0.5)
)

SweepLambda <- function(n, Lambda)
  colSums(sweep(Lambda, 2, n, "*"))

## -----------------------------------------------------------------------------

## *** prepare obs and sim data ***
samples <- c("V4", "V8") #row.names(moles.obs)

obs <- as.matrix(moles.obs[samples,-(1:2)])

sim <- matrix(0, length(samples), ncol(moles.sim), , list(samples, names(moles.sim)))
for(.sample in samples) {
  sim[.sample,] <- colMeans(moles.sim[index(moles.sim) >= moles.obs[.sample, "start"] & index(moles.sim) < moles.obs[.sample, "end"],])
}

## -----------------------------------------------------------------------------

## *** obs ***

j <- colnames(obs)
atomr.obs <- list()
omocm1.obs <- list()
for(.x in names(lambdaC)) {
  nC <- obs %*% lambdaC[[.x]][j]
  atomr.obs[[.x]] <- obs %*% t(Lambda[,j]) / nC[,1]
  omocm1.obs[[.x]] <- t(apply(obs, 1, SweepLambda, am[rownames(Lambda)]*Lambda[,j])) / (nC[,1]*am["C"])
}

## -----------------------------------------------------------------------------

## *** sim ***

cmpds <- intersect(colnames(sim), rownames(Y))
fgsets <- list(set1=intersect(meas.set1, colnames(Theta)), full=colnames(Theta))

atomr.sim <- list()
omocm1.sim <- list()
for(.x in names(fgsets)) {
  j <- fgsets[[.x]]
  k <- names(which(apply(Theta[,j] > 0, 1, any)))
  ctype.sim <- sim[,cmpds] %*% Y[cmpds,k]
  atomr.sim[[.x]] <- ctype.sim %*% Theta[k,j] %*% t(Lambda[,j]) / rowSums(ctype.sim)
  tmp <- t(apply(ctype.sim %*% Theta[k,j], 1, SweepLambda, am[rownames(Lambda)] *Lambda[,j])) / (rowSums(ctype.sim)*am["C"])
  tmp[,"CO"] <- tmp[,"CO"]+tmp[,"aldehyde"]
  omocm1.sim[[.x]] <- tmp[,-grep("aldehyde", colnames(tmp))]
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## *** van Krevlin diagram ***

pdf("outputs/fig_Sax.pdf", width=9.5, height=4.5)

colors.sample <- colorRampPalette(c("darkblue", "cornflowerblue"))(length(samples))
                                        #c("orange", "midnightblue")

par(mfrow=c(1, 2))
par(mar=c(2, 4, 2, 1), oma=c(2.5, 0, 0, 0))
par(mgp=c(1.8, .2, 0), tck=0.025, las=1)
par(cex=1.2)

pt.cex <- 1.2
plot.new()
plot.window(xlim=c(0, 1.2), ylim=c(0, 2.05), xaxs="i", yaxs="i")
for(.slope in c(0, -.5, -1, -2))
abline(2, .slope, lty=2)
points(atomr.obs[["old"]][,"O"], atomr.obs[["old"]][,"H"], col=colors.sample, cex=pt.cex, lwd=1.5)
points(atomr.obs[["new"]][,"O"], atomr.obs[["new"]][,"H"], pch=19, col=colors.sample, cex=pt.cex)
## arrows(atomr.obs[["old"]][,"O"], atomr.obs[["old"]][,"H"], #col=colors.sample,
##        atomr.obs[["new"]][,"O"], atomr.obs[["new"]][,"H"],
##        pch=19, code=2, angle=45, length=0.15)
points(atomr.sim[["set1"]][,"O"], atomr.sim[["set1"]][,"H"], pch=17, col=colors.sample, cex=pt.cex)
points(atomr.sim[["full"]][,"O"], atomr.sim[["full"]][,"H"], pch=15, col=colors.sample, cex=pt.cex)
axis(1)
axis(2)
box()
legend("bottomleft", pch=c(1, 19, 17, 15), legend=c("Obs., base case", "Obs., revised", "Sim., Set1", "Sim., Full"), bty="n", cex=.95)
legend(.45, par("usr")[3], xjust=0, yjust=0, legend=c("8h", "20h"), border=NA, fill=colors.sample, bty="n")
title(xlab="O/C", ylab="H/C", xpd=NA, cex.lab=1.2)
text(par("usr")[1], par("usr")[4]+par("cxy")[2]*.3, adj=c(0, 0), xpd=NA, "a)")

## *** OM/OC plot ***

group.levs <- c("aCH", "aCOH", "COOH", "CO", "CONO2", "eCH", "hydroperoxide", "peroxyacyl nitrate")

Stackrect <- function(x, y, ...) {
  ## lexically scoped: groups.levs, colors.FG, dx
  j <- intersect(group.levs, names(y))
  col <- colors.FG[j]
  dx <- dx
  ##
  positions <<- c(positions, x)
  ##
  yc <- c(0, cumsum(y[j]))
  for(i in seq(length(yc)-1))
    rect(x-dx, yc[i], x+dx, yc[i+1], col=col[i], ...)
}

dx <- .45
plot.new()
plot.window(c(0, 10), c(0, 1.4), yaxs="i")
positions <- c()
Stackrect(1, omocm1.obs[["old"]]["V4",],  border=NA)
Stackrect(2, omocm1.obs[["new"]]["V4",],  border=NA)
Stackrect(3, omocm1.sim[["set1"]]["V4",], border=NA)
Stackrect(4, omocm1.sim[["full"]]["V4",], border=NA)
Stackrect(6, omocm1.obs[["old"]]["V8",],  border=NA)
Stackrect(7, omocm1.obs[["new"]]["V8",],  border=NA)
Stackrect(8, omocm1.sim[["set1"]]["V8",], border=NA)
Stackrect(9, omocm1.sim[["full"]]["V8",], border=NA)
axis(1, positions, FALSE)
text(positions, par("usr")[3]-par("cxy")[2]*.5, adj=c(1, .5), srt=30,
     rep(c("Obs., base case", "Obs., revised", "Sim., Set1", "Sim., Full"), 2),
     xpd=NA)
yval <- seq(0, par("usr")[4], .2)
axis(2, yval, sprintf("%.1f", yval+1))
box()
text(c(2.5, 7.5), 2.3-1, adj=c(.5, .5), c(expression(underline("8h")), expression(underline("20h"))), xpd=NA)
title(ylab="OM/OC", cex.lab=1.2)
text(par("usr")[1], par("usr")[4]+par("cxy")[2]*.3, adj=c(0, 0), xpd=NA, "b)")

dev.off()

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
