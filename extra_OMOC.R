
options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(reshape2)
library(tidyverse)
library(Rfunctools)
library(RJSONIO)
library(gridExtra)
library(libRST)
theme_set(theme_bw())


PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_carbonprod.R", "lib/lib_metrics.R"))

carbon.attr <- ReadFile("carbonattr")
molec.attr <- ReadFile("molecattr")
## molec.ppb <- ReadFile("tseries_aer")
DBind[X, Y, Theta, gamma, zFG, Lambda] <- c(ReadFile("matrices"), ReadFile("matrices_2"))

molec.micromolm3 <- ReadMicromolm3(FilePath("tseries_aer"))

## -----------------------------------------------------------------------------

mw <- with(molec.attr, setNames(MW, compound))

## -----------------------------------------------------------------------------

melt.zoo <- function(data, varnames = names(dimnames(data)), ..., na.rm = FALSE, value.name = "value", index.include = TRUE, index.name = if(index.include) "index.time" else NULL) {
  mat <- as.matrix(data)
  .index <- setNames(index(data), rownames(mat))
  lf <- melt(mat, varnames = varnames, ..., na.rm = na.rm, as.is = TRUE, value.name = value.name)
  if(index.include) {
    .index.time <- setNames(sprintf("T%d", seq(nrow(mat))), rownames(mat))
    lf[[index.name]] <- .index.time[lf[[varnames[1]]]]
  }
  lf[[varnames[1]]] <- .index[lf[[varnames[1]]]]
  lf[c(index.name, varnames, value.name)]
}

## -----------------------------------------------------------------------------

conc.lf <- melt(molec.micromolm3[,rownames(X)], c("hour", "compound"), value.name = "micromolm3")

x.lf <- melt(X, c("compound", "fg"), as.is=TRUE, value.name="x") %>%
  filter(x > 0)

y.lf <- data.frame(compound = rownames(Y), nC = rowSums(Y), row.names=NULL)

lambda.lf <- melt(Lambda, c("element", "fg"), as.is = TRUE, value.name = "lambda") %>%
  filter(element != "C")

fg.lf <- inner_join(
  conc.lf,
  inner_join(x.lf, lambda.lf, by="fg"),
  by="compound"
)

carbon.lf <- inner_join(
  conc.lf, y.lf, by="compound"
)

sum.lf <- full_join(
  ## non-carbon atoms
  fg.lf %>%
  mutate(
    mole = micromolm3 * x * lambda,
    mass = mole * am[element]
  ) %>%
  group_by(index.time, hour, fg) %>%
  summarize(
    Oj = sum(mole[element=="O"]),
    OMj.m1 = sum(mass)
  ) %>%
  ungroup,
  ## carbon atoms
  carbon.lf %>%
  group_by(index.time, hour) %>%
  summarize(C = sum(micromolm3 * nC)) %>%
  mutate(OC = C * am["C"]) %>%
  ungroup,
  ##
  by = c("index.time", "hour")
) %>%
  mutate(OMOC.m1 = OMj.m1/OC, O.C = Oj/C)

## -----------------------------------------------------------------------------

if(FALSE) {

  ## check ratio of OM/OC to O/C: This is correct.

  sum.lf %>% group_by(fg) %>%
    mutate(ratj = OMj.m1 / Oj) %>%
    summarize(ratj = median(ratj, na.rm=TRUE)) %>%
    data.frame

  (t(Lambda[-1,]) %*% am[rownames(Lambda)[-1]]) / Lambda["O",]

}

## -----------------------------------------------------------------------------

group.levs <- rev(c("aCH", "aCOH", "COOH", "ketone", "aldehyde", "CONO2", "eCH", "hydroperoxide", "peroxyacyl nitrate"))

forplot <- sum.lf %>% filter(OMOC.m1 > 0)

forplot$fglabel <- Relabel(forplot$fg, labels.FG)

forplot <- forplot %>% filter(fglabel %in% group.levs) %>%
  mutate(fglabel = factor(fglabel, group.levs))

ggp <- list()

ggp[["O/C"]] <- ggplot(forplot) +
  geom_area(aes(hour, O.C, group=fglabel, fill=fglabel))+
  scale_fill_manual(name="Functional group", values = colors.FG[group.levs])+
  scale_x_continuous(name="Hour", expand=c(0, 0))+
  scale_y_continuous(name="O/C", limits = c(0, .8), expand=c(0, 0))+
  guides(fill = guide_legend(reverse=TRUE))

ggp[["OM/OC"]] <- ggplot(forplot) +
  geom_area(aes(hour, OMOC.m1, group=fglabel, fill=fglabel))+
  scale_fill_manual(name="Functional group", values = colors.FG[group.levs])+
  scale_x_continuous(name="Hour", expand=c(0, 0))+
  scale_y_continuous(name="OM/OC",
                     limits=c(1, 2.4)-1, breaks=seq(1, 2.4, .2)-1, labels=seq(1, 2.4, .2),
                     expand=c(0, 0)) +
  guides(fill = guide_legend(reverse=TRUE))


pdf("extra_outputs/extra_OMOC.pdf", width=6, height=5)
grid.arrange(GGMirrorTicks(ggp[["O/C"]]),
             GGMirrorTicks(ggp[["OM/OC"]]),
             nrow=2)
dev.off()

## -----------------------------------------------------------------------------

hour <- index(molec.micromolm3)

nC <- molec.micromolm3[,rownames(Y)] %*% as.matrix(rowSums(Y))
nC.2 <- molec.micromolm3[,rownames(X)] %*% X %*% as.matrix(Lambda["C",colnames(X)])

plot(hour, nC, type="l")
lines(hour, nC.2, type="l")

## -----------------------------------------------------------------------------

OM <- molec.micromolm3[,names(mw)] %*% mw
OM.2 <- sum.lf %>% group_by(index.time, hour) %>%
  summarize(OM = sum(OMj.m1) + OC[1]) %>% ungroup %>%
  arrange(hour)

## OM.3 <- molec.micromolm3[,rownames(X)] %*% X %*% sweep(t(Lambda), 2, am[rownames(Lambda)], "*")

plot(hour, OM, type="l", ylim=c(0, 1000))
with(OM.2, lines(hour, OM, col=2))
## lines(hour, rowSums(OM.3), col=3)

## -----------------------------------------------------------------------------
