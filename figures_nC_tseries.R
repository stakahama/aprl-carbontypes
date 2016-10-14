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
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_units.R", "lib/lib_metrics.R", "lib/lib_collapse.R"))
GGTheme()

## -----------------------------------------------------------------------------

matrices <- c(ReadFile("matrices"), ReadFile("matrices_2"))
molec.attr <- ReadFile("molecattr")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
clabels <- ReadFile("clabels")
carbon.attr <- ReadFile("carbonattr")
molec.moles <- ReadMicromolm3(FilePath("tseries_aer"))
decisions <- as.list(ReadFile("example_1"))
## lambdaC <- ReadFile("lambdaC_coef_actual")

ff <- list.files("outputs", "lambdaC_values.+\\.csv$", full=TRUE)
names(ff) <- sub(".+(set[0-9])(_collapsed)?\\.csv", "\\1\\2", basename(ff))

lambdaC <- ReadFile("lambdaC_array")

## -----------------------------------------------------------------------------

uniq <- UniqueMapping(matrices$Theta)

## -----------------------------------------------------------------------------

time <- index(molec.moles)
n <- coredata(molec.moles)

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

margins <- c("meas", "method")
nC.h <- do.call(
  Map, c(list(NCratio),
         do.call(expand.grid, c(dimnames(lambdaC)[margins], stringsAsFactors=FALSE)),
         list(MoreArgs=NamedList(lambdaC, n, matrices, measlist, time)))
)

tables <- ldply(unname(nC.h))

## -----------------------------------------------------------------------------

## tables <- ldply(tables, .id="case") %>%
##   melt(., c("time", "case"), variable.name="meas")

ggp <- ggplot(tables)+
  geom_hline(yintercept=1, size=.3, linetype=2)+
  geom_line(aes(time, value, color=meas))+
  facet_grid(.~method)+
  labs(x="Hour", y="Ratio")+
  lims(y=1+.2*c(-1,1))

## pdf(FilePath("plot_nC_tseries"))
## print(ggp)
## dev.off()

methods <- c("count", "solve", "fit", "nominal")
tables$method <- factor(tables$method, methods, labels.method[methods]) #toupper(methods))
tables$meas <- with(tables, factor(meas, unique(meas), Capitalize(unique(meas))))

## another viable option:
##
## ggplot(merged)+
##   geom_hline(yintercept=1, size=.2, color="gray")+
##   geom_line(aes(time, value, color=method, linetype=method))+
##   facet_wrap(~meas)+
##   lims(y=1+.2*c(-1,1))+
##   labs(y=expression(hat(italic(n))[C]^"*"/ italic(n)[C]^"*"))

lett.labels <- with(tables, {
  data.frame(letter=sprintf("%s) %s", letters[seq(nlevels(method))], levels(method)),
             method=factor(levels(method), levels(method)))
})

ex <- .2
ggp <- ggplot(tables)+
  ## geom_hline(yintercept=seq(1-ex, 1+ex, .05), linetype=2, size=.2, color="gray")+
  geom_hline(yintercept=1, size=.2, color="gray")+
  geom_line(aes(time, value, color=meas, linetype=meas), size=.8)+
  facet_wrap(~method)+
  lims(y=1+ex*c(-1,1))+
  labs(x="Hour", y=expression("Ratio,"~hat(italic(n))[C]^"*"/ italic(n)[C]^"*"))+
  scale_color_discrete(guide=guide_legend(title=""))+
  scale_linetype_discrete(guide=guide_legend(title=""))+
  geom_text(aes(x=-Inf, y=Inf, label=letter), data=lett.labels, size=5, hjust=0, vjust=1.2)+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

pdf(FilePath("plot_nC_recovery"), width=6, height=5)
print(ggp)
dev.off()

## tables %>% group_by(method, meas) %>% do(data.frame(value=round(with(., range(value))-1, 2)))
