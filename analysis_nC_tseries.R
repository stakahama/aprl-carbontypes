options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_io.R", "lib/lib_metrics.R"))

## -----------------------------------------------------------------------------

matrices <- c(ReadFile("matrices"), ReadFile("matrices_2"))

molec.attr <- ReadFile("molecattr")

DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]

clabels <- ReadFile("clabels")

carbon.attr <- ReadFile("carbonattr")

moles.molec <- ReadTSeries(FilePath("tseries_aer"))

decisions <- as.list(ReadFile("example_1"))

lambdaC <- list(nominal=ReadFile("lambdaC_coef_nominal"),
                actual=ReadFile("lambdaC_coef_actual"))

DBind[X, Y, Theta, gamma, zFG, Lambda] <- matrices

## -----------------------------------------------------------------------------

cmpds <- intersect(names(moles.molec), rownames(Y))

time <- index(moles.molec)
n <- coredata(moles.molec)

nC <- n[,cmpds] %*% rowSums(Y[cmpds,])
fg <- n[,cmpds] %*% X[cmpds,]

## -----------------------------------------------------------------------------

tables <- list()

for(.est in names(lambdaC)) {

  lambdaC[[.est]] <- replace(lambdaC[[.est]], is.na(lambdaC[[.est]]), 0)

  nC.h <- list()

  for(.label in names(measlist)) {

    j <- measlist[[.label]]

    Y.s <- sweep(Y, 2, sign(rowSums(Theta[,j]*gamma[j])), `*`)
    nC.s <- n[,cmpds] %*% rowSums(Y.s[cmpds,])

    nC.h[[.label]] <- (fg[,j] %*% lambdaC[[.est]][.label,j]) / nC.s

  }

  tables[[.est]] <- cbind(time=time, as.data.frame(nC.h))

}

## -----------------------------------------------------------------------------

tables <- ldply(tables, .id="case") %>%
  melt(., c("time", "case"), variable.name="meas")

tables$case <- factor(tables$case, c("actual", "nominal"))

ggp <- ggplot(tables)+
  geom_hline(yintercept=1, size=.3, linetype=2)+
  geom_line(aes(time, value, color=meas))+
  facet_grid(case~.)+
  labs(x="Hour", y="Ratio")+
  lims(y=1+.2*c(-1,1))

pdf(FilePath("plot_nC_tseries"))
print(ggp)
dev.off()
