

options(stringsAsFactors=FALSE)

library(abind)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(xtable)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", "lib/lib_lambdaC.R")

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")
lambdaCdf <- ReadFile("lambdaC_coef_actual_f")
decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

ReplZero <- function(x) replace(x, is.na(x), 0)

Nominal <- function (lambdaC, Theta, adjust = NULL) {
    nom <- c(1/sort(unique(rowSums(Theta))), 0)
    lambdaC[] <- sapply(lambdaC, function(x, y) if(is.na(x)) x else y[which.min(abs(x -y))], nom)
    if (!is.null(adjust))
        for (j in names(adjust)) lambdaC[, j] <- adjust[j]
    lambdaC
}


id <- c("method", "meas")
method.order <- c("count", "solve", "fit")
lamClist <- dlply(lambdaCdf, .(method), function(x) acast(x, meas~group, value.var="Estimate"))
lamCnom <- Nominal(lamClist[["count"]], Theta, decisions$lambdaC.nominal)
lamCarr <- abind(c(lamClist[method.order], list(nominal=lamCnom)), along=length(dim(lamCnom))+1)
names(dimnames(lamCarr)) <- c(id[2], "group", id[1])
lambdaC <- aperm(lamCarr, c(1,3,2))

saveRDS(lambdaC, FilePath("lambdaC_array"))

## -----------------------------------------------------------------------------

method.order <- dimnames(lambdaC)$method
fixed.lookup <- with(unique(subset(lambdaCdf, ,c(group, fixed))), setNames(fixed, group))

dfm <- full_join(lambdaCdf,
                 melt(lambdaC[,"nominal",], c("meas", "group"), value.name="Estimate") %>%
                 mutate(method="nominal", fixed=fixed.lookup[as.character(group)]))

dfm$text <- with(dfm, {
  ifelse(is.na(Estimate), "",
  ifelse(is.na(`Std. Error`), sprintf("%.2f", Estimate),
         sprintf("%.2f (%.2f)", Estimate, `Std. Error`)))
})

wf <- dcast(dfm %>% filter(!fixed) %>% select(meas, method, group, text),
            meas+method~group, value.var="text", fill="") %>%
  mutate(meas=factor(meas), method=factor(method, method.order)) %>%
  arrange(meas, method)

levels(wf$meas) <- Capitalize(levels(wf$meas))
levels(wf$method) <- toupper(levels(wf$method))

print(xtable(wf[,c("meas", "method", "alkane CH", "alcohol", "organonitrate", "alkene CH", "hydroperoxide")]),
      include.rownames=FALSE, sanitize.text.function = identity)
