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

## -----------------------------------------------------------------------------

ReplZero <- function(x) replace(x, is.na(x), 0)

margins <- c("meas", "method")
lambdaC <- daply(lambdaC, margins,
                 function(x, j) ReplZero(unlist(x[setdiff(names(x), j)])),
                 margins)

## -----------------------------------------------------------------------------

time <- index(molec.moles)
n <- coredata(molec.moles)

## nC <- n[,cmpds] %*% rowSums(Y[cmpds,])
## fg <- n[,cmpds] %*% X[cmpds,]

## -----------------------------------------------------------------------------


NCratio <- function(meas, method, lambdaC, n, matrices, measlist, time) {
  DBind[X, Y, Theta, gamma] <- matrices
  ##
  i <- intersect(colnames(n), rownames(Y))
  j <- intersect(measlist[[meas]], colnames(Theta))
  ##
  G <- n[,i] %*% X[i,]
  Y.s <- sweep(Y, 2, sign(Theta[,j] %*% gamma[j]), `*`)
  nC.s <- n[,i] %*% rowSums(Y.s[i,])
  nC.r <- (G[,j] %*% lambdaC[meas, method, j]) / nC.s
  ##
  data.frame(time, meas, method, value=nC.r)
}

nC.h <- do.call(
  Map, c(list(NCratio),
         do.call(expand.grid, c(dimnames(lambdaC)[margins], stringsAsFactors=FALSE)),
         list(MoreArgs=NamedList(lambdaC, n, matrices, measlist, time)))
)

tables <- ldply(unname(nC.h))

## tables <- ldply(tables, .id="case") %>%
##   melt(., c("time", "case"), variable.name="meas")

ggp <- ggplot(tables)+
  geom_hline(yintercept=1, size=.3, linetype=2)+
  geom_line(aes(time, value, color=meas))+
  facet_grid(.~method)+
  labs(x="Hour", y="Ratio")+
  lims(y=1+.2*c(-1,1))

pdf(FilePath("plot_nC_tseries"))
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

## *** collapsed ***

DBind[X, Y, Theta, gamma] <- matrices

cmpds <- intersect(colnames(n), rownames(Y))

ave.nom <- Nominal(lambdaC, "count", decisions$lambdaC.nominal)

NCratio2 <- function(meas, lambdaC, n, X, Y, Theta, gamma, measlist, time) {
  j <- intersect(measlist[[meas]], colnames(Theta))
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
                                 measlist, time)))
)

tables2 <- ldply(unname(nC.h))

ggp <- ggplot(tables2)+
  geom_line(aes(time, value, color=meas))+
  facet_wrap(~method)

## print(ggp)

## -----------------------------------------------------------------------------

merged <- full_join(
  tables,
  tables2 %>% mutate(method="nominal")
)


methods <- c("count", "solve", "fit", "nominal")
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

lett.labels <- with(merged, data.frame(letter=sprintf("%s) %s", letters[seq(nlevels(method))], levels(method)),
                               method=factor(levels(method), levels(method))))

ggp <- ggplot(merged)+
  geom_hline(yintercept=1, size=.2, color="gray")+
  geom_line(aes(time, value, color=meas, linetype=meas))+
  facet_wrap(~method)+
  lims(y=1+.2*c(-1,1))+
  labs(x="Hour", y=expression("Recovery fraction,"~hat(italic(n))[C]^"*"/ italic(n)[C]^"*"))+
  scale_color_discrete(guide=guide_legend(title=""))+
  scale_linetype_discrete(guide=guide_legend(title=""))+
  geom_text(aes(x=-Inf, y=Inf, label=letter), data=lett.labels, size=5, hjust=0, vjust=1.2)+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

pdf("outputs/nC_recovery_tseries.pdf", width=6, height=5)
GGTheme()
print(ggp)
dev.off()

