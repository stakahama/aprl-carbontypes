options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_metrics.R"))

## -----------------------------------------------------------------------------

matrices <- c(ReadFile("matrices"), ReadFile("matrices_2"))
molec.attr <- ReadFile("molecattr")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
clabels <- ReadFile("clabels")
carbon.attr <- ReadFile("carbonattr")
molec.moles <- ReadMicromolm3(FilePath("tseries_aer"))
decisions <- as.list(ReadFile("example_1"))
lambdaC <- list(actual=ReadFile("lambdaC_coef_actual_f"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma, zFG, Lambda] <- matrices

zeta <- with(unique(carbon.attr[,c("ctype", "OSC")]), setNames(OSC, ctype))[rownames(Theta)]
cOM <- Ctypemass(Theta, gamma, Lambda)

n.example <- coredata(Slice(molec.moles, decisions$hour))[1,]
cmpds <- intersect(names(n.example), rownames(Y))
carbons <- rownames(Theta)

groups <- setdiff(colnames(X), grep("radical", colnames(X), value=TRUE))

full <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda, zeta, cOM, cmpds, carbons, groups)

amsOSC.est <- with(full$ratios, 2*sum(`O/C`) - sum(`H/C`))

## -----------------------------------------------------------------------------

Multijoin <- function(x, y1, y2, elem) {
  x <- x[[elem]] %>% mutate(meas="full", case="ideal")
  y1 <- ldply(y1, "[[", elem, .id="meas") %>% mutate(case="ideal")
  y2 <- ldply(y2, "[[", elem, .id="meas") %>% mutate(case="estimated")
  full_join(full_join(x, y1), y2) %>%
    mutate(case=factor(case, c("ideal", "estimated")))
}

CumsumFrac <- function(x)
  cumsum(sort(x, decreasing=TRUE))/sum(x)

PlotCumsums <- function(mat) {
  par(mfrow=n2mfrow(ncol(mat)), oma=c(1, 1, 0, 0))
  for(v in colnames(mat)) {
    cumsum. <- CumsumFrac(mat[,v])
    xval <- seq_along(cumsum.)
    plot(xval, cumsum., ylim=c(0, 1), xaxt="n", type="o",
         main=v, xlab="", ylab="cumulative sum")
    axis(1, xval, FALSE)
    text(xval, par("usr")[3]-par("cxy")[2]*.5, names(cumsum.), adj=c(1, .5), srt=30, xpd=NA)
  }
}

CumsumDF <- function(x, index)
  ResetIndex(data.frame(cumsum=CumsumFrac(x)), index)

for(.est in names(lambdaC)) {

  out <- out.h <- list()
  ## .h is the version where carbon number is estimated

  for(.label in names(measlist)) {

    meas <- intersect(measlist[[.label]], colnames(Theta))

    ## ctypes <- intersect(rownames(Theta)[rowSums(Theta[,meas]) > 0], carbons)
    ctypes <- rownames(Theta)[rowSums(Theta[,meas]) > 0]

    ## star <- Calculate(n.example, X, Y, Theta, gamma, zFG, Lambda, cmpds, ctypes, meas)
    star <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda,
                       zeta, cOM, cmpds, ctypes, meas)
    star.h <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda,
                         zeta, cOM, cmpds, ctypes, meas,
                         lambdaC[[.est]][.label,])

    out[[.label]] <- star
    out.h[[.label]] <- star.h

  }

  ## -----------------------------------------------------------------------------

  merged.c <- Multijoin(full, out, out.h, "masses") %>%
    mutate(clabel=c(clabels, setNames(,"X"))[ctype], ctype=NULL) %>%
    melt(., c("clabel", "meas", "case")) %>%
    group_by(variable) %>%
    arrange(clabel)

  merged.g <- Multijoin(full, out, out.h, "ratios") %>%
    melt(., c("group", "meas", "case")) %>%
    group_by(variable) %>%
    arrange(group)

  merged.osc <- Multijoin(full, out, out.h, "OSC") %>%
    full_join(., data.frame(method="approx", meas="AMS", value=amsOSC.est, case=factor("ideal", c("ideal", "estimated")))) %>%
    mutate(method=factor(method, c("true", "approx")))

  saveRDS(list(mass=merged.c, props=merged.g, OSC=merged.osc), SprintF("props_file", .est))

}
