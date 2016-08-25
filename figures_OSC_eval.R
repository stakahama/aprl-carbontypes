
options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R", "lib/lib_units.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
DBind[,,merged.osc] <- ReadFile(SprintF("props_file", "actual"))
molec.moles <- lapply(setNames(FilePath(c("tseries_gas", "tseries_aer")), c("gas", "aer")),
                      ReadMicromolm3)
clabels <- ReadFile("clabels")
carbon.attr <- ReadFile("carbonattr")
molec.attr <- ReadFile("molecattr")
decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

examp <- ldply(molec.moles, function(X, i, Y)
  CarbonProd(Slice(X, i), Y) %>% mutate(time=NULL),
  decisions$hour, Y, .id="phase")

## -----------------------------------------------------------------------------

ctypes <- ldply(measlist, function(j, Theta)
  data.frame(ctype=rownames(Theta)[rowSums(Theta[,intersect(j, colnames(Theta))]) > 0]),
  Theta, .id="meas")

wf <- full_join(left_join(ctypes, examp %>% select(phase, ctype, nC)),
                examp %>% select(phase, ctype, nC) %>% mutate(meas="full")) %>%
  dcast(., phase+ctype~meas, sum, value.var="nC")

meas. <- c(names(measlist), "full")
wf[,meas.] <- cbind(wf[,meas.[1],drop=FALSE], t(apply(wf[,meas.], 1, diff)))

## -----------------------------------------------------------------------------

nC <- left_join(melt(wf, measure.vars=meas., variable.name="meas", value.name="nC"),
                carbon.attr %>% distinct(ctype, OSC))

osc <- merged.osc %>%
  ungroup() %>% filter(case=="ideal") %>% mutate(case=NULL) %>%
  filter((meas=="full" & method=="true") | (meas!="full" & method=="approx")) %>%
  mutate(method=NULL)

## -----------------------------------------------------------------------------

Addbars <- function(height, at, col, ...) {
  x <- 0
  for(i in rownames(height)) {
    x2 <- x+height[i,]
    rect(x, at[,1], x2, at[,2], col=col[i], border=NA, ...)
    x <- x2
 }
}

levels(nC$meas) <- Capitalize(levels(nC$meas))

mat <- acast(nC %>% filter(phase=="aer"), meas~OSC, sum, value.var="nC")
colors.sets <- setNames(brewer.pal(8, "Set1"), rownames(mat))
osc$index <- seq_along(osc$meas)
osc$meas <- Capitalize(osc$meas)
yval <- as.integer(colnames(mat))

xmax <- ceiling(max(colSums(mat)))

pdf("outputs/OSC_distr_meas.pdf", width=8, height=4)
## layout(t(1:2), width=c(3,2))
layout(t(1:2))
par(tck=0.025, mgp=c(1.8, .2, 0), oma=c(2,2,1,1))
par(mar=c(2,2,.5,.5))
## par(mar=c(2,2,.5,4))
ylim <- c(-3, 3)
i <- 1
##
plot.new()
plot.window(c(0, xmax), ExpandLim(ylim), xaxs="i")
ypos <- outer(seq(min(yval), max(yval)), .35*c(-1,1), "+") # fixed
Addbars(mat, at=ypos, col=colors.sets)
axis(1)
axis(2, las=1)
box()
## legend(par("usr")[2], par("usr")[4], xjust=0, yjust=1, xpd=NA,
##        legend=rownames(mat), fill=col, border=NA, bty="n")
legend(par("usr")[2], par("usr")[4], xjust=1, yjust=1, xpd=NA,
       legend=rownames(mat), fill=colors.sets[rownames(mat)], border=NA, bty="n")
mtext(expression(italic(n)[C]~(mu*"mole"/m^3)), 1, line=par("mgp")[1])
mtext(expression(OS[C]), 3)
text(par("usr")[1], par("usr")[4], adj=c(0, -.3), xpd=NA,
     sprintf("%s)", letters[i]), cex=1.2)
i <- i+1
##
par(mar=c(2,2,.5,.5))
with(osc, {
  plot.new()
  plot.window(ExpandLim(range(index), .1), ExpandLim(ylim))
  abline(h=seq(min(ylim), max(ylim)), lty=2, col=8)
  abline(h=0)
  lines(index, value, type="h", lwd=2, col="midnightblue")
  points(index, value, pch=19, lwd=2, cex=1.2, col="midnightblue")
  axis(1, index, FALSE)
  axis(2, las=1)
  axis(3,,FALSE)
  axis(4,,FALSE)
  box()
  text(index, par("usr")[3]-par("cxy")[2]*.3, adj=c(1, .5), xpd=NA, srt=30,
       ifelse(meas=="AMS", expression(2*O/C-H/C), meas))
  mtext(expression(bar(OS)[C]*"*"), 3)
})
text(par("usr")[1], par("usr")[4], adj=c(0, -.3), xpd=NA,
     sprintf("%s)", letters[i]), cex=1.2)
mtext("Carbon oxidation state", 2, outer=TRUE)#las=0, line=par("mgp")[1], outer=TRUE)
dev.off()
