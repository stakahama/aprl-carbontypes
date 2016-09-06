
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
carbon.attr <- ReadFile("carbonattr")

## -----------------------------------------------------------------------------

Lambda <- Lambda[c("H", "O", "N"),]
colnames(Lambda) <- Relabel(colnames(Lambda), c(labels.FG, setNames("CO", "ketone"))) # We will use the column of ketone for CO as it dominates the condensed phase over aldehyde (the only difference is in the no. of H). For formal aggregation, use AggGroups() in lib/lib_collapse.R, which assumes equal molar portions the mean.
colnames(Theta) <- colnames(Lambda)
meas.set1 <- Relabel(measlist$set1, c(labels.FG, setNames("CO", "ketone")))
names(zFG) <- Relabel(names(zFG), c(labels.FG, setNames("CO", "ketone")))

lambdaC <- list(
  old=c("aCH"=.63, "COOH"=1, "CO"=1, "aCOH"=0.63, "CONO2"=0.63),
  new=c("aCH"=.45, "COOH"=1, "CO"=1, "aCOH"=0.5, "CONO2"=0.5)
)

SweepLambda <- function(n, Lambda)
  colSums(sweep(Lambda, 2, n, "*"))

## -----------------------------------------------------------------------------

## *** prepare obs and sim data ***
samples <- c("V2", "V8") #row.names(moles.obs)

obs <- as.matrix(moles.obs[samples,-(1:2)])

sim <- matrix(0, length(samples), ncol(moles.sim), , list(samples, names(moles.sim)))
for(.sample in samples) {
  sim[.sample,] <- colMeans(moles.sim[index(moles.sim) >= moles.obs[.sample, "start"] & index(moles.sim) < moles.obs[.sample, "end"],])
}

## par(mfrow=c(1, 2))
## pie(obs["V2",])
## pie(obs["V8",])

## -----------------------------------------------------------------------------

## *** obs ***

j <- colnames(obs)
atomr.obs <- list()
omocm1.obs <- list()
osc.obs <- list()
for(.x in names(lambdaC)) {
  nC <- (obs %*% lambdaC[[.x]][j])[,1]
  atomr.obs[[.x]] <- obs %*% t(Lambda[,j]) / nC
  omocm1.obs[[.x]] <- t(apply(obs, 1, SweepLambda, am[rownames(Lambda)]*Lambda[,j])) / (nC*am["C"])
  osc.obs[[.x]] <- (obs %*% zFG[j])[,1]/nC
}

## -----------------------------------------------------------------------------

## *** sim ***

cmpds <- intersect(colnames(sim), rownames(Y))
fgsets <- list(set1=intersect(meas.set1, colnames(Theta)), full=colnames(Theta))
zeta <- with(unique(subset(carbon.attr,,c(ctype, OSC))), setNames(OSC, ctype))

atomr.sim <- list()
omocm1.sim <- list()
osc.sim <- list()
for(.x in names(fgsets)) {
  j <- fgsets[[.x]]
  k <- names(which(apply(Theta[,j] > 0, 1, any)))
  ctype.sim <- sim[,cmpds] %*% Y[cmpds,k]
  nC <- rowSums(ctype.sim)
  atomr.sim[[.x]] <- ctype.sim %*% Theta[k,j] %*% t(Lambda[,j]) / nC
  tmp <- t(apply(ctype.sim %*% Theta[k,j], 1, SweepLambda, am[rownames(Lambda)] *Lambda[,j])) / (nC*am["C"])
  tmp[,"CO"] <- tmp[,"CO"]+tmp[,"aldehyde"]
  omocm1.sim[[.x]] <- tmp[,-grep("aldehyde", colnames(tmp))]
  osc.sim[[.x]] <- (ctype.sim %*% zeta[k])[,1]/nC
}

atomr.apin <- local({
  nC <- sum(Y["APINENE",])
  (Y["APINENE",] %*% Theta %*% t(Lambda))[1,]/nC
})

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

labels.est <- c("MEAS-PREV", "MEAS-NOM", "SIM-SET1", "SIM-FULL")

pdf(FilePath("plot_Sax"), width=14, height=5)

par(mfrow=c(1, 3))
par(mar=c(2, 4, 2, 1), oma=c(2.5, 0, 0, 0))
par(mgp=c(1.8, .2, 0), tck=0.025, las=1)
par(cex=1.2)

colors.sample <- colorRampPalette(c("darkblue", "cornflowerblue"))(length(samples))
                                        #c("orange", "midnightblue")

## *** van Krevlin diagram ***

pt.cex <- 1.2
plot.new()
plot.window(xlim=c(0, 1.2), ylim=c(0, 2))#, xaxs="i", yaxs="i")
for(.slope in c(0, -1, -2))
  abline(2, .slope, lty=2)
points(atomr.apin["O"], atomr.apin["H"], pch=18, cex=1.6, xpd=NA)
points(atomr.sim[["set1"]][,"O"], atomr.sim[["set1"]][,"H"], pch=15, col=colors.sample, cex=pt.cex)
points(atomr.sim[["full"]][,"O"], atomr.sim[["full"]][,"H"], pch=17, col=colors.sample, cex=pt.cex)
points(atomr.obs[["old"]][,"O"], atomr.obs[["old"]][,"H"], col=colors.sample, cex=pt.cex, lwd=1.5)
points(atomr.obs[["new"]][,"O"], atomr.obs[["new"]][,"H"], pch=21, col="white", bg=colors.sample, cex=pt.cex)
## arrows(atomr.obs[["old"]][,"O"], atomr.obs[["old"]][,"H"], #col=colors.sample,
##        atomr.obs[["new"]][,"O"], atomr.obs[["new"]][,"H"],
##        pch=19, code=2, angle=45, length=0.15)
axis(1)
axis(2)
box()
legend("bottomleft", pch=c(18, 1, 19, 15, 17), pt.cex=c(1.4, 1, 1, 1, 1),
       legend=c(expression(alpha*"-pinene"), labels.est),
       bty="n", cex=.95)
legend(.5, par("usr")[3], xjust=0, yjust=0, legend=c("4h", "21h"), border=NA, fill=colors.sample, bty="n")
title(xlab="O/C", ylab="H/C", xpd=NA, cex.lab=1.2)
text(par("usr")[1], par("usr")[4]+par("cxy")[2]*.3, adj=c(0, 0), xpd=NA, "a)", cex=1.2)

## *** OM/OC plot ***

group.levs <- c("aCH", "aCOH", "COOH", "CO", "CONO2", "eCH", "hydroperoxide", "peroxyacyl nitrate")

Stackrect <- function(x, y, ...) {
  ## lexically scoped: groups.levs, colors.FG, dx
  ## global var: positions
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

early <- samples[1]
late <- samples[2]
dx <- .45
plot.new()
plot.window(c(0, 10), c(0, 1.4), yaxs="i")
positions <- c()
Stackrect(1, omocm1.obs[["old"]][early,],  border=NA)
Stackrect(2, omocm1.obs[["new"]][early,],  border=NA)
Stackrect(3, omocm1.sim[["set1"]][early,], border=NA)
Stackrect(4, omocm1.sim[["full"]][early,], border=NA)
Stackrect(6, omocm1.obs[["old"]][late,],   border=NA)
Stackrect(7, omocm1.obs[["new"]][late,],   border=NA)
Stackrect(8, omocm1.sim[["set1"]][late,],  border=NA)
Stackrect(9, omocm1.sim[["full"]][late,],  border=NA)
axis(1, positions, FALSE)
text(positions, par("usr")[3]-par("cxy")[2]*.5, adj=c(1, .5), srt=30,
     rep(labels.est, 2),
     xpd=NA)
yval <- seq(0, par("usr")[4], .2)
axis(2, yval, sprintf("%.1f", yval+1))
box()
text(c(2.5, 7.5), par("usr")[4]-par("cxy")[2]*.5, adj=c(.5, 1), c(expression(underline("4h")), expression(underline("21h"))), xpd=NA, cex=1.1)
title(ylab="OM/OC", cex.lab=1.2)
text(par("usr")[1], par("usr")[4]+par("cxy")[2]*.3, adj=c(0, 0), xpd=NA, "b)", cex=1.2)


htype <- function(x, y, ...) {
  positions <<- c(positions, x)
  lines(x, y, type="h", lwd=2, col="midnightblue", ...)
  points(x, y, pch=19, cex=1.2, lwd=2, col="midnightblue", ...)
}

early <- samples[1]
late <- samples[2]
plot.new()
plot.window(c(0, 10), c(-3, 3))
abline(h=seq(-3, 3), lty=2, col=8)
abline(h=0)
positions <- c()
htype(1, osc.obs[["old"]][early])
htype(2, osc.obs[["new"]][early])
htype(3, osc.sim[["set1"]][early])
htype(4, osc.sim[["full"]][early])
htype(6, osc.obs[["old"]][late])
htype(7, osc.obs[["new"]][late])
htype(8, osc.sim[["set1"]][late])
htype(9, osc.sim[["full"]][late])
axis(1, positions, FALSE)
text(positions, par("usr")[3]-par("cxy")[2]*.5, adj=c(1, .5), srt=30,
     rep(labels.est, 2),
     xpd=NA)
axis(2)
box()
text(c(2.5, 7.5), par("usr")[4]-par("cxy")[2]*.5, adj=c(.5, 1), c(expression(underline("4h")), expression(underline("21h"))), xpd=NA, cex=1.1)
title(ylab=expression(bar(OS)[C]))
text(par("usr")[1], par("usr")[4]+par("cxy")[2]*.3, adj=c(0, 0), xpd=NA, "c)", cex=1.2)

dev.off()

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
