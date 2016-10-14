
################################################################################
##
## figures_props.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################


options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)
library(ggplot2)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_C_attributes.R", "lib/lib_OSC.R", "lib/lib_metrics.R"))
## source("http://ms.mcmaster.ca/~bolker/R/misc/legendx.R")
source("lib/legendx.R")

## -----------------------------------------------------------------------------

DBind[merged.c, merged.g, merged.osc] <- ReadFile(SprintF("props_file", "actual"))

## -----------------------------------------------------------------------------

SelectCase <- function(x)
  x %>% ungroup() %>% filter(case=="ideal") %>% mutate(case=NULL)

## -----------------------------------------------------------------------------

## mass

mass <- SelectCase(merged.c)

ref <- with(mass %>% filter(meas=="full") %>% group_by(variable) %>% summarize(value=sum(value)),
            setNames(value, variable))

mass$value[] <- with(mass, value/ref[as.character(variable)])

cumul <- mass %>%
  filter(meas=="full" & variable=="OC") %>% mutate(meas=NULL, variable=NULL) %>%
  arrange(desc(value)) %>% mutate(value=round(cumsum(value), 2))

## with(cumul, plot(value, type="o"))

clabel.levs <- with(cumul, clabel[1:20])

mass <- mass %>% filter(clabel %in% clabel.levs) %>%
  mutate(clabel=factor(clabel, sort(as.integer(clabel.levs))))

## ggplot(mass)+
##   geom_bar(aes(meas, value, fill=clabel), stat="identity")+
##   facet_grid(.~variable)

## -----------------------------------------------------------------------------

df <- SelectCase(merged.g)

atomr <- df %>% filter(variable!="OM/OC-1")
omoc <- df %>% filter(variable=="OM/OC-1")

cumul <- omoc %>%
  filter(meas=="full") %>% mutate(meas=NULL, variable=NULL) %>%
  arrange(desc(value)) %>% mutate(value=round(cumsum(value)/sum(value), 2))

## group.keep <- with(cumul, group[1:8])

group.levs <- c("alkane CH", "alcohol", "carboxylic acid", "ketone", "aldehyde", "organonitrate", "aromatic CH", "alkene CH", "phenol", "hydroperoxide", "peroxyacylnitrate")

## call factor twice to drop unused levels
atomr <- atomr %>% filter(group %in% group.levs) %>%
  mutate(group=factor(factor(group, group.levs)))
omoc <- omoc %>% filter(group %in% group.levs) %>%
  mutate(group=factor(factor(group, group.levs)))

## -----------------------------------------------------------------------------


## meas.levs <- with(mass, setNames(sort(unique(meas)),  Capitalize(sort(unique(meas)))))
## mass$meas <- factor(mass$meas, meas.levs, names(meas.levs))
## atomr$meas <- factor(atomr$meas, meas.levs, names(meas.levs))
## omoc$meas <- factor(omoc$meas, meas.levs, names(meas.levs))

mass$meas <- factor(mass$meas, names(labels.meas), labels.meas)
atomr$meas <- factor(atomr$meas, names(labels.meas), labels.meas)
omoc$meas <- factor(omoc$meas, names(labels.meas), labels.meas)

## -----------------------------------------------------------------------------

LegendFig <- function(...) {
  plot.new()
  plot.window(0:1, 0:1)
  get("legend", globalenv())(..., bty="n", xpd=NA)
}

barplot <- function(...) {
  dotargs <- list(...)
  dotargs$xaxt <- "n"
  dotargs$border <- NA
  bx <- do.call(graphics::barplot, c(dotargs, list(yaxt="n")))
  axis(1,bx,FALSE)
  if(is.null(dotargs$yaxt))
    axis(2, cex.axis=1.4)
  axis(3,bx,FALSE)
  axis(4,,FALSE)
  box()
  text(bx, par("usr")[3]-par("cxy")[2]*.3, adj=c(1, .5),
       colnames(dotargs[[1]]), xpd=NA, srt=30, cex=1.4)
  text(par("usr")[1], par("usr")[4],
       sprintf("%s)", letters[i]), xpd=NA, adj=c(0, -.3), cex=1.4)
  i <<- i+1
}

Parset <- function() {
  par(mar=c(2, 2, 2, 1), oma=c(1, 2, 0, 0))
  par(mgp=c(1.8, .2, 0), tck=0.025, las=1)
}

ylims <- list(
  "OC"=c(0, 1.02),
  "OM"=c(0, 1.02),
  "O/C"=c(0, 1.02),
  "H/C"=c(0, 2.02),
  "N/C"=c(0, 0.05),
  "OM/OC"=c(0, 1.02)
)

colors.C <- with(list(x=clabel.levs), setNames(GGColorHue(length(x)), x))
cex.exp <- c(1.2, .8)

pdf(FilePath("plot_props_1"), width=7, height=2.9)
layout(matrix(1:3, ncol=3, byrow=TRUE), width=c(1, 1, .5))
Parset()
i <- 1
##
for(.var in levels(mass$variable)) {
  .table <- filter(mass, variable==.var)
  .mat <- acast(.table, clabel~meas, fill=0)
  barplot(.mat[names(colors.C),], ylim=ylims[[.var]], col=colors.C)
  mtext(.var, 3)
}
mtext("Recovery fraction", 2, outer=TRUE, las=0, line=.5)
with(list(x=names(colors.C)),
     LegendFig(x=.35, xjust=.5, y=.5, yjust=.5, title="Carbon type", legend=x, fill=colors.C[x], ncol=2, border=NA, box.cex=cex.exp))
dev.off()

pdf(FilePath("plot_props_2"), width=10.5, height=2.5)
layout(matrix(1:5, ncol=5, byrow=TRUE), c(rep(1, 4), .7))
Parset()
i <- 1
##
for(.var in c("O/C", "H/C", "N/C")) {
  .table <- filter(atomr, variable==.var)
  .mat <- acast(.table, group~meas, fill=0)
  barplot(.mat, ylim=ylims[[.var]], col=colors.FG[Relabel(rownames(.mat),labels.FG)])
  mtext(.var, 3)
}
##
.mat <- acast(omoc, group~meas, fill=0)
barplot(.mat, ylim=ylims[["OM/OC"]], yaxt="n", col=colors.FG[Relabel(rownames(.mat), labels.FG)])
mtext("OM/OC", 3)
yval <- seq(0, par("usr")[4], .2)
axis(2, yval, sprintf("%.1f", yval+1), cex.axis=1.4)
## mtext(c("Ratio", "Recovery fraction"), 2, adj=c(.24, .8), outer=TRUE, las=0)
mtext("Ratio", 2, outer=TRUE, las=0, line=.5)
with(list(x=Relabel(levels(omoc$group),labels.FG)),
     LegendFig(x=.1, xjust=0.5, y=.5, yjust=.5, title="FG", legend=x, fill=colors.FG[x], ncol=1, border=NA, box.cex=cex.exp, text.width=.3))
dev.off()

## -----------------------------------------------------------------------------

osc <- SelectCase(merged.osc) %>%
   filter((meas=="full" & method=="true") | (meas!="full" & method=="approx")) %>%
   mutate(method=NULL)

extra <- setdiff(osc$meas, names(labels.meas))
osc$meas <- factor(osc$meas, c(names(labels.meas), extra), c(labels.meas, extra))
osc$index <- unclass(osc$meas)#seq(nrow(osc))

pdf(FilePath("plot_osc_f"), width=7, height=5)
par(mfrow=c(1,1), cex=1.2)
Parset()
with(osc, {
  plot.new()
  plot.window(range(index), c(-4, 3), yaxs="i")
  abline(h=seq(-4, 3), lty=2, col=8)
  abline(h=0)
  lines(index, value, type="h", lwd=2, col="midnightblue")
  points(index, value, pch=19, lwd=2, col="midnightblue")
  axis(1, index, FALSE)
  axis(2)
  axis(3,,FALSE)
  axis(4,,FALSE)
  box()
  text(index, par("usr")[3]-par("cxy")[2]*.3, adj=c(1, .5), xpd=NA, srt=30,
       ifelse(meas=="AMS", expression(2*O/C-H/C), as.character(meas)))
  mtext(expression(bar(OS)[C]), 2, las=0, line=par("mgp")[1])
})
dev.off()
