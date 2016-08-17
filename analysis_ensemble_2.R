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

n.example <- coredata(Slice(moles.molec, decisions$hour))[1,]
cmpds <- intersect(names(which(n.example > 1e-3)), rownames(Y))

groups <- setdiff(names(which(apply(X[cmpds,] > 0, 2, any))),
                  grep("radical", colnames(X), value=TRUE))

## full <- Calculate(n.example, X, Y, Theta, gamma, zFG, Lambda, cmpds, , groups)
full <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda, carbon.attr, cmpds, , groups)

## -----------------------------------------------------------------------------

Multijoin <- function(x, y1, y2, elem) {
  x <- x[[elem]] %>% mutate(meas="full", case="ideal")
  y1 <- ldply(y1, "[[", elem, .id="meas") %>% mutate(case="ideal")
  y2 <- ldply(y2, "[[", elem, .id="meas") %>% mutate(case="estimated")
  full_join(full_join(x, y1), y2) %>%
    mutate(case=factor(case, c("ideal", "estimated")))
}

for(.est in names(lambdaC)) {

  out <- out.h <- list()
  ## .h is the version where carbon number is estimated

  for(.label in names(measlist)) {

    meas <- measlist[[.label]]

    ctypes <- rownames(Theta)[rowSums(Theta[,meas]) > 0]
    ## star <- Calculate(n.example, X, Y, Theta, gamma, zFG, Lambda, cmpds, ctypes, meas)
    star <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda,
                       carbon.attr, cmpds, ctypes, meas)
    star.h <- Calculate2(n.example, X, Y, Theta, gamma, zFG, Lambda,
                         carbon.attr, cmpds, ctypes, meas, lambdaC[[.est]][.label,])

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

  saveRDS(list(mass=merged.c, props=merged.c), SprintF(FilePath("props_file"), .est))

  wf <- acast(merged.c %>% filter(case=="ideal"), variable~meas, sum)
  recovery <- wf/wf[,"full"]

  print(recovery)

  wf.c <- dcast(merged.c, variable+meas~case, sum) %>%
    filter(estimated > 0)

  wf.g <- dcast(merged.g, variable+meas~case, sum) %>%
    filter(estimated > 0)

  ## -----------------------------------------------------------------------------

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

  ggp.c <- ggplot(merged.c)+
    facet_grid(variable~case)+
    geom_bar(aes(meas, value, fill=clabel), stat="identity")

  ggp.g <- ggplot(merged.g)+
    facet_grid(variable~case, scale="free_y")+
    geom_bar(aes(meas, value, fill=group), stat="identity")

  pdf(SprintF("plot_props", .est), width=16, height=7)
  print(grid.arrange(ggp.c, ggp.g, ncol=2))
  dev.off()

}
