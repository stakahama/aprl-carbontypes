

options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(plotrix)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R", "lib/lib_units.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma, zFG, Lambda] <- c(ReadFile("matrices"), ReadFile("matrices_2"))
molec.moles <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadMicromolm3)
clabels <- ReadFile("clabels")
examp.FG <- ReadFile("fgorder")
carbon.attr <- ReadFile("carbonattr")
molec.attr <- ReadFile("molecattr")
decisions <- as.list(ReadFile("example_1"))
measlist <- ReadFile("meas")$lambdaC

## -----------------------------------------------------------------------------

RevDict <- function(x)
  setNames(rep(names(x), sapply(x, length)), unlist(x))

RecDiff <- function(x)
  setNames(
    Reduce(function(init, x)
      c(init, list(setdiff(x, unlist(init)))),
      x, c()),
    names(x)
  )

diffs.FG <- RecDiff(c(measlist, list(full=colnames(Theta))))
names(diffs.FG)[] <- Capitalize(names(diffs.FG))

sets.ctype <- lapply(measlist, function(j, Theta)
  names(which(apply(Theta[,intersect(j, colnames(Theta))]>0, 1, any))),
  Theta)
diffs.ctype <- RecDiff(c(sets.ctype, list(full=rownames(Theta))))
names(diffs.ctype)[] <- Capitalize(names(diffs.ctype))

diffs.FG <- RevDict(diffs.FG)
diffs.ctype <- RevDict(diffs.ctype)

names(diffs.FG)[] <- Relabel(names(diffs.FG), labels.FG)
names(diffs.ctype)[] <- clabels[names(diffs.ctype)]

## -----------------------------------------------------------------------------

first <- function(var, x, margin) {
  if(missing(margin)) {
    var==head(names(x), 1)
  } else {
    var==head(dimnames(x)[[margin]], 1)
  }
}

last <- function(var, x, margin) {
  if(missing(margin)) {
    var==tail(names(x), 1)
  } else {
    var==tail(dimnames(x)[[margin]], 1)
  }
}

## -----------------------------------------------------------------------------

Zoolf <- function(x, varnames=c("index", "name")) {
  xval <- setNames(index(x), rownames(x))
  df <- melt(coredata(x), varnames=varnames, as.is=TRUE)
  df[[varnames[1]]] <- xval[df[[varnames[1]]]]
  df
}

sel <- lapply(molec.moles, Slice, decisions$hour)
sel <- ldply(sel, Zoolf, varnames=c("index", "compound"), .id="phase")

## cmpds <- intersect(sel$compound, rownames(Y))

ctype <- inner_join(sel %>% mutate(index=NULL),
                    melt(Y, c("compound", "ctype"), value.name="count", as.is=TRUE) %>% mutate(nC=rowSums(Y)[compound]),
                    by="compound") %>%
  group_by(phase, nC, ctype) %>%
  summarize(value=sum(count*value))

ctype <- ctype %>%
  group_by(phase) %>%
  mutate(value=value/sum(value))

ctype <- left_join(ctype, unique(subset(carbon.attr,,c(ctype, OSC))))
ctype$clabel <- clabels[as.character(ctype$ctype)]

grid <- unclass(by(ctype, ctype[c("phase", "OSC")], function(x)
  acast(x, nC~clabel)))

grid <- grid[,ncol(grid):1]

## -----------------------------------------------------------------------------

mycol <- colorRampPalette(c("white", "darkred"))(64)

zlim <- c(0, .22)
zval <- seq(zlim[1], zlim[2], ,length(mycol))
zlab.val <- with(list(z=pretty(zval)), z[findInterval(z, range(zlim))==1])
zlab <- rep("", length(zval))
iz <- sapply(zlab.val, function(x, y) which.min(abs(y-x)), zval)
zlab[iz] <- sprintf("%.2f", zlab.val)

pdf(FilePath("plot_nCgrid_ctype"), width=7, height=7)
nC <- unique(ctype$nC)
layout(cbind(matrix(seq(length(grid)), c(ncol(grid), nrow(grid))), length(grid)+1),
       width=c(5, 5, 1), height=sapply(grid["gas",], ncol))
par(oma=c(4, 10, 2, 2))
par(mar=c(.2,.2,.2,.2), mgp=c(2, .2, 0), tck=0.02)
cex.axis <- 1.5
for(.phase in rownames(grid)) {
  for(.osc in colnames(grid)) {
    mat <- grid[[.phase, .osc]]
    xval <- seq_along(nC)
    yval <- seq(ncol(mat))
    image(xval, yval, mat[rev(as.character(nC)),rev(yval),drop=FALSE], axes=FALSE, ann=FALSE,
          col=mycol, zlim=zlim)
    line <- 5
    dx <- 0
    if(first(.phase, grid, 1)) {
      axis(2, yval, FALSE) #rev(colnames(mat)), las=1, cex.axis=cex.axis)
      text(par("usr")[1]-par("cxy")[1]*.3, yval, rev(colnames(mat)), adj=c(1, .5),
           cex=cex.axis, col=colors.set[diffs.ctype[rev(colnames(mat))]], xpd=NA)
      if(first(.osc, grid, 2)) {
        mtext(bquote(OS[C]==.(.osc)), 2, las=1, line=line, xpd=NA)
      } else {
        mtext(bquote(.(.osc)), 2, las=1, line=line, xpd=NA)
      }
      arrows(par("usr")[1]-par("usr")[1]*(line-dx), par("usr")[3],
             par("usr")[1]-par("usr")[1]*(line-dx), par("usr")[4],
             code=3, angle=90, length=0.01, xpd=NA)
    }
    if(first(.osc, grid, 2)) {
      mtext(c("gas"="Gas", "aer"="Aerosol")[.phase], 3, line=.5)
    }
    if(last(.osc, grid, 2)) {
      ## axis(1, xval, rev(nC), xpd=NA, cex.axis=cex.axis)
      text(xval, par("usr")[3]-par("cxy")[2]*.1, adj=c(.5, 1), rev(nC),
           cex=cex.axis, xpd=NA)
    }
    for(.side in c(1, 3)) axis(.side, xval, FALSE)
    for(.side in c(2, 4)) axis(.side, yval, FALSE)
    box()
  }
}
mtext(expression(italic(n)[C]), 1, adj=.45, outer=TRUE, line=1.8)
mtext("Carbon type", 2, outer=TRUE, line=1.8)
plot.new()
plot.window(0:1, 0:1)
xleft <- .2
dx <- .2
color.legend(xleft, 0, xleft+dx, 1, zlab, rect.col=mycol, cex=.8, align="rb", gradient="y")
mtext("Carbon type fraction in each phase", 4, outer=TRUE, line=.5)
dev.off()

## -----------------------------------------------------------------------------

Tmp <- function(x, lab, OSC) {
  ix <- rownames(x)
  rownames(x) <- lab[ix]
  split(as.data.frame(x), OSC[ix])
}

thet <- lapply(Tmp(Theta, clabels, with(unique(subset(carbon.attr,,c(ctype, OSC))), setNames(OSC, ctype))), function(x) { x <- t(as.matrix(x)); x[,order(as.numeric(colnames(x))),drop=FALSE] })

thet <- rev(thet)

fg <- with(list(x=zFG[colnames(Theta)], n=examp.FG[colnames(Theta)]),
           names(x[order(-x, n)]))#Relabel(names(x), labels.FG))]))
xval <- seq_along(fg)
xpos <- rev(sapply(split(xval, zFG[fg]), mean))
xrg <- sapply(split(xval, zFG[fg]), range)
xrg <- xrg[,ncol(xrg):1]

## mycol <- colorRampPalette(c("white", "darkviolet"))(max(Theta)+1)
mycol <- colorRampPalette(c("white", "black"))(max(Theta)+1)
zlim <- range(Theta)

pdf(FilePath("plot_nCgrid_fg"), width=7, height=7)
layout(cbind(matrix(seq(length(thet))), length(thet)+1),
       width=c(5, 1), height=sapply(thet, ncol))
par(oma=c(1, 10, 14, 1))
par(mar=c(.2,.2,.2,.2), mgp=c(2, .2, 0), tck=0.02)
cex.axis <- 1.5
for(.osc in names(thet)) {
  mat <- thet[[.osc]]
  yval <- seq(ncol(mat))
  image(xval, yval, mat[fg,rev(yval),drop=FALSE], axes=FALSE, ann=FALSE,
        col=mycol, zlim=zlim)
  axis(2, yval, FALSE) #rev(colnames(mat)), las=1, cex.axis=cex.axis)
  text(par("usr")[1]-par("cxy")[1]*.3, yval, rev(colnames(mat)), adj=c(1, .5),
       cex=cex.axis, col=colors.set[diffs.ctype[rev(colnames(mat))]], xpd=NA)
  ## ----------
  line <- 6
  dx <- 1
  if(first(.osc, thet)) {
    mtext(bquote(OS[C]==.(.osc)), 2, xpd=NA, line=line-dx, las=1)
  } else {
    mtext(.osc, 2, xpd=NA, line=line-dx, las=1)
  }
  arrows(par("usr")[1]-par("usr")[1]*line, par("usr")[3],
         par("usr")[1]-par("usr")[1]*line, par("usr")[4],
         code=3, angle=90, length=0.01, xpd=NA)
  ## ----------
  line <- .3
  expf <- 8.7
  dx <- 0.0
  if(first(.osc, thet)) {
    fg.newlab <- Relabel(fg, labels.FG)
    text(xval, par("usr")[4]+par("cxy")[2]*line, adj=c(0, .5), srt=60, xpd=NA,
         fg.newlab, cex=cex.axis, col=colors.set[diffs.FG[fg.newlab]])
    ## segments(xrg[1,]-dx, par("usr")[4]+par("cxy")[2]*(expf-.5),
    ##          xrg[2,]+dx, par("usr")[4]+par("cxy")[2]*(expf-.5), xpd=NA)
    arrows(xrg[1,]-dx, par("usr")[4]+par("cxy")[2]*(expf-.5),
           xrg[2,]+dx, par("usr")[4]+par("cxy")[2]*(expf-.5),
           code=3, angle=90, length=0.01, xpd=NA)
    text(xpos, par("usr")[4]+par("cxy")[2]*expf,
         c(bquote(italic(z)==.(names(xpos)[1])), parse(text=names(xpos)[-1])),
         xpd=NA, cex=cex.axis)
  }
  ## ----------
  for(.side in c(1, 3)) axis(.side, xval, FALSE)
  axis(4, yval, FALSE)
  box()
}
mtext("Carbon type", 2, outer=TRUE, line=2.5)
##
plot.new()
plot.window(0:1, 0:1)
xleft <- .2
dx <- .2
color.legend(xleft, 0, xleft+dx, 1, seq_along(mycol)-1, rect.col=mycol, cex=.8, align="rb", gradient="y")
mtext("Number of groups associated with each carbon type", 4, outer=TRUE, line=-2)
dev.off()
