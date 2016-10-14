
################################################################################
##
## figures_nC_unc.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## see LICENSE_GPLv3.txt
##
################################################################################


options(stringsAsFactors=FALSE)

library(Rfunctools)
library(RJSONIO)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
## PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R", "lib/lib_units.R"))
GGTheme()

## -----------------------------------------------------------------------------

OMOC <- seq(1, 2.5, .1)
OCr <- seq(0, 1.2, .1)
unc <- .1

## -----------------------------------------------------------------------------

f <- function(OMOC, delta) 1-1/(1+delta)-1/OMOC+1/(OMOC*(1+delta))
f2 <- function(OCr, delta) 1-1/(1+delta)

if(FALSE) {

  par(mfrow=c(1, 2))

  plot.new()
  plot.window(c(1, 2.5), c(1, 2.5))
  axis(1)
  axis(2)
  box()
  matlines(OMOC, 1+outer(OMOC-1, 1+unc*c(-1,1), "/"), type="l", lty=2)

  plot.new()
  plot.window(c(0, 1.2), c(0, 1.2))
  axis(1)
  axis(2)
  box()
  matlines(OCr, outer(OCr, (1+unc*c(-1,1)), "/"), type="l", lty=2)

  ## -----------------------------------------------------------------------------

  par(mfrow=c(2, 2))

  delta <- c(0.05, .1)
  matplot(OMOC, 1-(1+outer(OMOC-1, 1+delta, "/"))/OMOC, type="l")

  matplot(OMOC, outer(OMOC, delta, f), type="l")

  matplot(OCr, 1-outer(OCr, 1+delta, "/")/OCr, type="l", ylim=c(0, .1))

  matplot(OCr, outer(OCr, delta, f2), type="l", ylim=c(0, .1))

}

## -----------------------------------------------------------------------------

delta <- seq(0.0, .1, .01)
names(delta) <- sprintf("X%d", seq(length(delta)))

wf <- rbind(data.frame(metric="OM/OC", nominal=OMOC, 1-(1+outer(OMOC-1, 1+delta, "/"))/OMOC),
            data.frame(metric='italic(n)[a]^"*"/italic(n)[C]^"*"', nominal=OCr, outer(OCr, delta, f2)))
lf <- melt(wf, c("metric", "nominal"), variable.name="delta", as.is=TRUE)
lf$delta <- delta[lf$delta]

limits <- rbind(data.frame(metric="OM/OC", x=c(1, 2.5), y=c(0, .1)),
                data.frame(metric='italic(n)[a]^"*"/italic(n)[C]^"*"', x=c(0, 1.201), y=c(0, .1)))

## *** attempt to make 'delta' symbol italic (does not export well to PDF) ***

## delta.Utf8 <- intToUtf8(letters.greek[4,"hexUTF"])

## ggp <- ggplot(lf)+
##   geom_blank(aes(x, y), data=limits)+
##   geom_line(aes(nominal, value, color=delta, group=delta))+
##   facet_grid(.~metric, scale="free_x", labeller = label_parsed)+
##   scale_x_continuous(expand=c(0, 0))+
##   scale_y_continuous(limits=c(0, .1))+
##   scale_color_continuous(name=bquote(italic(.(delta.Utf8))["["*italic(n)[C]^"*"*"]"]), limits=c(0, .1))+
##   labs(x="Nominal value", y=bquote(italic(.(delta.Utf8))))+
##   theme(panel.margin.x = unit(1.5, "lines"))

## png("outputs/nC_unc_deltas.png", width=7, height=3.5, units="in", res=96)
## print(ggp)
## dev.off()

## *** use built-in 'delta' symbol ***

panels <- data.frame(metric=c('italic(n)[a]^"*"/italic(n)[C]^"*"', "OM/OC"), label=c("a)", "b)"))

ggp <- ggplot(lf)+
  geom_blank(aes(x, y), data=limits)+
  geom_text(aes(-Inf, Inf, label=label), hjust=-.3, vjust=1.2, data=panels, size=5)+
  geom_line(aes(nominal, value, color=delta, group=delta))+
  facet_grid(.~metric, scale="free_x", labeller = label_parsed)+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(limits=c(0, .1))+
  scale_color_continuous(name=expression(delta["["*italic(n)[C]^"*"*"]"]), limits=c(0, .1))+
  labs(x="Actual value", y=expression("Relative error,"~delta))+
  theme(panel.margin.x = unit(1.5, "lines"))

pdf(FilePath("plot_nC_unc"), width=7, height=3.5)
print(ggp)
dev.off()
