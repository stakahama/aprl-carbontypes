
################################################################################
##
## analysis_IMPROVE.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################


options(stringsAsFactors=FALSE)

library(zoo)
library(chron)
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
DBind[, measlist, collapserule] <- ReadFile("meas")
carbon.attr <- ReadFile("carbonattr")

path <- "../../../Multivariate/plsdev_runs/fitspath_FGOMLV"

.sites <- c("PMRFX", "OLYMX", "SAMAX")
.months <- list(MAM=c("Mar", "Apr", "May"), JJA=c("Jun", "Jul", "Aug"))

pred <- read.delim(file.path(path, "run_9101/M2c_amb_prediction_table.txt"))
carbox <- read.delim(file.path(path, "run_9102/M2c_amb_prediction_table.txt"))
pred <- full_join(pred %>% filter(!variable %in% c("CO", "cCOH")), carbox)
pred[,c("site", "date")] <- do.call(rbind, strsplit(pred$sample, "_"))
pred$date <- as.chron(pred$date, "%Y%m%d")

## -----------------------------------------------------------------------------

sitedate <- SetIndex(unique(pred[,c("sample", "site", "date")]), "sample")
obs <- acast(pred %>% filter(site %in% .sites & months(date) %in% unlist(.months)),
             sample~variable, value.var="predicted")


## -----------------------------------------------------------------------------

Lambda <- Lambda[c("H", "O", "N"),]
colnames(Lambda) <- Relabel(colnames(Lambda), c(labels.FG, setNames("CO", "ketone"))) # We will use the column of ketone for CO as it dominates the condensed phase over aldehyde (the only difference is in the no. of H). For formal aggregation, use AggGroups() in lib/lib_collapse.R, which assumes equal molar portions the mean.
colnames(Theta) <- colnames(Lambda)
meas.set1 <- Relabel(measlist$set1, c(labels.FG, setNames("CO", "ketone")))
names(zFG) <- Relabel(names(zFG), c(labels.FG, setNames("CO", "ketone")))

lambdaC <- list(
  old=c("aCH"=.5, "COOH"=1, "CO"=1, "aCOH"=0),
  new=c("aCH"=.45, "COOH"=1, "CO"=1, "aCOH"=0.5)
)

## -----------------------------------------------------------------------------

SweepLambda <- function(n, Lambda)
  colSums(sweep(Lambda, 2, n, "*"))

BoundNeg <- function(x) {
  x[] <- replace(x, x < 0, 0)
  x
}

j <- colnames(obs)

nC.obs <- list()
omocm1.obs <- list()
for(.x in names(lambdaC)) {
  nC <- (obs %*% lambdaC[[.x]][j])[,1]
  omocm1.obs[[.x]] <- t(apply(obs, 1, SweepLambda, am[rownames(Lambda)]*Lambda[,j])) / (nC*am["C"])
  omocm1.obs[[.x]] <- BoundNeg(omocm1.obs[[.x]])
  nC.obs[[.x]] <- nC
}

omocm1.table <- ldply(omocm1.obs, function(x, id)
  melt(cbind(id[rownames(x), ], x), c("site", "date"), variable.name="group"),
  sitedate, .id="version")

nC.table <- ldply(nC.obs, function(x)
  data.frame(sample=names(x), value=x), .id="version")

## -----------------------------------------------------------------------------

RevDict <- function(x)
  setNames(rep(names(x), times=sapply(x, length)), unlist(x))

omocm1.agg <- omocm1.table %>%
  mutate(season=RevDict(.months)[as.character(months(date))]) %>%
  group_by(version, site, season, group) %>%
  summarize(value=mean(value))

ggplot(omocm1.agg)+
  geom_bar(aes(site, value, fill=group), stat="identity")+
  facet_grid(season~version)

nC.wf <- dcast(nC.table, sample~version)
with(nC.wf, plot(old, new))
abline(0, 1)
