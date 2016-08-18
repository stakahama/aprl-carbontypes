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

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma, zFG, Lambda] <- matrices

zeta <- with(unique(carbon.attr[,c("ctype", "OSC")]), setNames(OSC, ctype))[rownames(Theta)]
cOM <- with(CarbonTypeMass(carbon.attr), setNames(OM, ctype))[rownames(Theta)]
amsOSC <- with(carbon.attr, setNames(2*O-H, ctype))[rownames(Theta)]

n.example <- coredata(Slice(moles.molec, decisions$hour))[1,]
cmpds <- intersect(names(n.example), rownames(Y))
## carbons <- names(sort(colSums(n.example[cmpds] * Y[cmpds,]), decreasing=TRUE)[decisions$ncarbons])
carbons <- rownames(Theta)

## groups <- setdiff(names(which(apply(X[cmpds,] > 0, 2, any))),
##                   grep("radical", colnames(X), value=TRUE))
## groups <- setdiff(names(which(apply(Theta[carbons,] > 0, 2, any))),
##                   grep("radical", colnames(X), value=TRUE))
groups <- setdiff(colnames(X), grep("radical", colnames(X), value=TRUE))

## full <- Calculate(n.example, X, Y, Theta, gamma, zFG, Lambda, cmpds, , groups)
full <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda, zeta, cOM, cmpds, carbons, groups)

amsOSC.est <- sum(Y[cmpds,] %*% amsOSC)/sum(Y[cmpds,])

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

    meas <- measlist[[.label]]

    ## ctypes <- intersect(rownames(Theta)[rowSums(Theta[,meas]) > 0], carbons)
    ctypes <- rownames(Theta)[rowSums(Theta[,meas]) > 0]

    ## star <- Calculate(n.example, X, Y, Theta, gamma, zFG, Lambda, cmpds, ctypes, meas)
    star <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda,
                       zeta, cOM, cmpds, ctypes, meas)
    star.h <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda,
                         zeta, cOM, cmpds, ctypes, meas, lambdaC[[.est]][.label,])

    out[[.label]] <- star
    out.h[[.label]] <- star.h

  }

  ## -----------------------------------------------------------------------------

  merged.c <- Multijoin(full, out, out.h, "masses") %>%
    mutate(clabel=clabels[ctype], ctype=NULL) %>%
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

  saveRDS(list(mass=merged.c, props=merged.g), SprintF("props_file", .est))

  ## -----------------------------------------------------------------------------

  agg.c <- acast(merged.c %>% filter(meas=="full"), clabel~variable, sum)

  agg.g <- acast(merged.g %>% filter(meas=="full"), group~variable, sum)

  saveRDS(list(mass=apply(agg.c, 2, CumsumDF, "clabel"), props=apply(agg.g, 2, CumsumDF, "group")),
          SprintF("tables_props_cumsum", .est))

  pdf(SprintF("plot_props_cumsum", .est), width=12, height=6)
  PlotCumsums(agg.c)
  PlotCumsums(agg.g)
  dev.off()

  ## -----------------------------------------------------------------------------

  ggp.c <- ggplot(merged.c)+
    facet_grid(variable~case)+
    geom_bar(aes(meas, value, fill=clabel), stat="identity")

  ggp.g <- ggplot(merged.g)+
    facet_grid(variable~case, scale="free_y")+
    geom_bar(aes(meas, value, fill=group), stat="identity")

  ggp.osc <- ggplot(merged.osc)+
    facet_grid(.~case)+
    geom_bar(aes(meas, value, fill=method), stat="identity", position="dodge")

  pdf(SprintF("plot_props", .est), width=20, height=7)
  print(grid.arrange(ggp.c, ggp.g, ggp.osc, ncol=3))
  dev.off()

  ## -----------------------------------------------------------------------------

  wf <- acast(merged.c %>% filter(case=="ideal"), variable~meas, sum)
  recovery <- wf/wf[,"full"]

  print(recovery)

  wf.c <- dcast(merged.c, variable+meas~case, sum) %>%
    filter(estimated > 0)

  wf.g <- dcast(merged.g, variable+meas~case, sum) %>%
    filter(estimated > 0)

  ggp.c <- ggplot(wf.c)+
    facet_wrap(~variable)+
    geom_point(aes(ideal, estimated, color=meas))+
    geom_abline(intercept=0, slope=1)

  ggp.g <- ggplot(wf.g)+
    facet_wrap(~variable)+
    geom_point(aes(ideal, estimated, color=meas))+
    geom_abline(intercept=0, slope=1)

  pdf(SprintF("plot_props_scatter", .est), width=16, height=7)
  print(grid.arrange(ggp.c, ggp.g, ncol=2))
  dev.off()

}
