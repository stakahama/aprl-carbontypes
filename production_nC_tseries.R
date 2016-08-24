options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(gridExtra)
library(ggplot2)
library(pryr)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_metrics.R", "lib/lib_collapse.R"))

## -----------------------------------------------------------------------------

matrices <- c(ReadFile("matrices"), ReadFile("matrices_2"))
molec.attr <- ReadFile("molecattr")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
clabels <- ReadFile("clabels")
carbon.attr <- ReadFile("carbonattr")
molec.moles <- ReadMicromolm3(FilePath("tseries_aer"))
decisions <- as.list(ReadFile("example_1"))
lambdaC <- ReadFile("lambdaC_coef_actual")

ff <- list.files("outputs", "lambdaC_values.+\\.csv$", full=TRUE)
names(ff) <- sub(".+(set[0-9])(_collapsed)?\\.csv", "\\1\\2", basename(ff))

## -----------------------------------------------------------------------------

uniq <- UniqueMapping(matrices$Theta)
uniq.collapsed <- UniqueMapping(matrices.collapsed$Theta)

for(j in intersect(colnames(lambdaC$actual), uniq$group)) {
  lambdaC$actual[,j] <- ifelse(is.na(lambdaC$actual[,j]), uniq[match(j, uniq$group),"value"], lambdaC$actual[,j])
  lambdaC$nominal[,j] <- ifelse(is.na(lambdaC$nominal[,j]), uniq[match(j, uniq$group),"value"], lambdaC$nominal[,j])
}

for(j in intersect(colnames(lambdaC$actual), uniq.collapsed$group)) {
  lambdaC$actual[,j] <- ifelse(is.na(lambdaC$actual[,j]), uniq[match(j, uniq$group),"value"], lambdaC$actual[,j])
  lambdaC$nominal[,j] <- ifelse(is.na(lambdaC$nominal[,j]), uniq[match(j, uniq$group),"value"], lambdaC$nominal[,j])
}

## -----------------------------------------------------------------------------

ReplZero <- function(x) replace(x, is.na(x), 0)

margins <- c("case", "meas", "method")
lambdaC <- daply(ldply(lambdaC, .id="case"), margins,
                 function(x, j) ReplZero(unlist(x[setdiff(names(x), j)])),
                 margins)

## -----------------------------------------------------------------------------

time <- index(molec.moles)
n <- coredata(molec.moles)

## nC <- n[,cmpds] %*% rowSums(Y[cmpds,])
## fg <- n[,cmpds] %*% X[cmpds,]

## -----------------------------------------------------------------------------


NCratio <- function(case, meas, method, lambdaC, n, arr1, arr2, measlist, time) {
  if(grepl("collapsed", meas)) {
    DBind[X, Y, Theta, gamma] <- arr2
  } else {
    DBind[X, Y, Theta, gamma] <- arr1
  }
  ##
  i <- intersect(colnames(n), rownames(Y))
  j <- measlist[[meas]]
  ##
  G <- n[,i] %*% X[i,]
  Y.s <- sweep(Y, 2, sign(Theta[,j] %*% gamma[j]), `*`)
  nC.s <- n[,i] %*% rowSums(Y.s[i,])
  nC.r <- (G[,j] %*% lambdaC[case, meas, method, j]) / nC.s
  ##
  data.frame(time, case, meas, method, value=nC.r)
}

nC.h <- do.call(
  Map, c(list(NCratio),
         do.call(expand.grid, c(dimnames(lambdaC)[1:3], stringsAsFactors=FALSE)),
         list(MoreArgs=NamedList(lambdaC, n,
                                 arr1=matrices, arr2=matrices.collapsed,
                                 measlist=c(measlist, measlist.collapsed),
                                 time)))
)

tables <- ldply(unname(nC.h))

## tables <- ldply(tables, .id="case") %>%
##   melt(., c("time", "case"), variable.name="meas")

tables$case <- factor(tables$case, c("actual", "nominal"))

ggp <- ggplot(tables)+
  geom_hline(yintercept=1, size=.3, linetype=2)+
  geom_line(aes(time, value, color=meas))+
  facet_grid(case~method)+
  labs(x="Hour", y="Ratio")+
  lims(y=1+.2*c(-1,1))

pdf(FilePath("plot_nC_tseries"))
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

## *** collapsed ***

DBind[X, Y, Theta, gamma] <- matrices.collapsed

cmpds <- intersect(colnames(n), rownames(Y))

nom <- c(1/sort(unique(rowSums(Theta))), 0)

ave <- apply(lambdaC["actual",grepl("collapsed", dimnames(lambdaC)[[2]]),,], c(1,3), mean)
ave.nom <- `[<-`(force(ave), , sapply(ave, function(x, y) y[which.min(abs(x-y))], nom))

NCratio2 <- function(meas, lambdaC, n, X, Y, Theta, gamma, measlist, time) {
  j <- measlist[[meas]]
  Y.s <- sweep(Y, 2, sign(Theta[,j] %*% gamma[j]), `*`)
  G <- n %*% X
  nC.s <- n %*% rowSums(Y.s)
  nC.r <- (G[,j] %*% lambdaC[meas, j]) / nC.s
  data.frame(time, meas, value=nC.r)
}

nC.h <- do.call(
  Map, c(list(NCratio2),
         dimnames(ave.nom)[1],
         list(MoreArgs=NamedList(lambdaC=ave.nom,
                                 n=n[,cmpds], X=X[cmpds,], Y=Y[cmpds,], Theta, gamma,
                                 measlist=measlist.collapsed, time)))
)

tables2 <- ldply(unname(nC.h))

ggp <- ggplot(tables2)+
  geom_line(aes(time, value, color=meas))+
  lims(y=1+.2*c(-1,1))

print(ggp)

## -----------------------------------------------------------------------------

merged <- full_join(
  tables %>% filter(grepl("\\_collapsed", meas) & case=="actual") %>%
  mutate(meas=sub("\\_collapsed", "", meas), case=NULL),
  tables2 %>% mutate(method="nominal", meas=sub("\\_collapsed", "", meas))
)


methods <- c("count", "solve", "nominal")
merged$method <- factor(merged$method, methods, toupper(methods))
merged$meas <- with(merged, factor(meas, unique(meas), Capitalize(unique(meas))))

## another viable option:
##
## ggplot(merged)+
##   geom_hline(yintercept=1, size=.2, color="gray")+
##   geom_line(aes(time, value, color=method, linetype=method))+
##   facet_wrap(~meas)+
##   lims(y=1+.2*c(-1,1))+
##   labs(y=expression(hat(italic(n))[C]^"*"/ italic(n)[C]^"*"))

lett.labels <- with(merged, data.frame(letter=sprintf("%s)", letters[seq(nlevels(method))]),
                                       method=factor(levels(method), levels(method))))

ggp <- ggplot(merged)+
  geom_hline(yintercept=1, size=.2, color="gray")+
  geom_line(aes(time, value, color=meas, linetype=meas))+
  facet_wrap(~method)+
  lims(y=1+.2*c(-1,1))+
  labs(x="Hour", y=expression("Recovery fraction,"~hat(italic(n))[C]^"*"/ italic(n)[C]^"*"))+
  scale_color_discrete(guide=guide_legend(title=""))+
  scale_linetype_discrete(guide=guide_legend(title=""))+
  geom_text(aes(x=-Inf, y=Inf, label=letter), data=lett.labels, size=5, hjust=0, vjust=1)

pdf("outputs/nC_recovery_tseries.pdf", width=8, height=3.5)
GGTheme()
print(ggp)
dev.off()

