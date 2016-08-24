
options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_collapse.R", "lib/lib_constrOptim.R", "lib/lib_lambdaC.R", "lib/lib_units.R"))

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadFile("matrices")
molec.attr <- ReadFile("molecattr")
molec.moles <- ReadMicromolm3(FilePath("tseries_aer"))
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]

## -----------------------------------------------------------------------------

cmpds <- intersect(names(molec.moles), rownames(Y))

## -----------------------------------------------------------------------------

Wrap <- function(y, X, fixed=c("carboxylic acid"=1, "carbonyl"=1, "ketone"=1, "phenol"=1, "peroxyacylnitrate"=1, "carbonylperoxyacid"=1)) {
  j.f <- intersect(colnames(X), names(fixed))
  j.inc <- setdiff(colnames(X), names(fixed))
  nC.fixed <- X[,j.f,drop=FALSE] %*% fixed[j.f]
  nC <- y - nC.fixed
  data <- list(y=nC, X=X[,j.inc,drop=FALSE])
  ##
  soln.lm <- lm(y~X-1, data=data)
  soln.blsq <- with(data, BoundedLsq(y, X))
  ##
  fn <- function(x) setNames(x, sub("^X","",names(x)))
  ## coef(summary(soln.lm))
  rbind(lm=c(fn(coef(soln.lm)), fixed[j.f])[colnames(X)],
        blsq=c(soln.blsq, fixed[j.f])[colnames(X)])
}

uniq <- with(UniqueMapping(Theta), setNames(value, group))

time <- index(molec.moles)
period <- TRUE# time > 2
fit <- list()
nCfit <- list()

for(.label in names(measlist)) {

  meas <- intersect(measlist[[.label]], colnames(Theta))
  meas <- meas[apply(X[cmpds,meas] > 0, 2, any)]

  measC <- rowSums(Theta[,meas]) > 0
  Y.s <- coredata(molec.moles)[period,cmpds] %*% Y[cmpds,measC]
  X.s <- coredata(molec.moles)[period,cmpds] %*% X[cmpds,meas]

  ## calculation of lambda

  ## fit.lm0 <- lm(y~X-1, data=list(y=Y.s, X=X.s))
  ## fit.lm <- lm(y~X-1, data=list(y=rowSums(Y.s), X=X.s))
  fit[[.label]] <- Wrap(y=rowSums(Y.s), X=X.s, uniq)
  soln <- fit[[.label]]["blsq",]
  nCfit[[.label]] <- data.frame(
    time=time,
    ref=rowSums(coredata(molec.moles)[,cmpds] %*% Y[cmpds,measC]),
    est=coredata(molec.moles)[,cmpds] %*% X[cmpds,meas] %*% soln
  )

}

cc <- ldply(fit, function(x) data.frame(as.list(x["blsq",]), check.names=FALSE), .id="meas")

nCfit <- ldply(nCfit, .id="meas")
nCfit$ratio <- with(nCfit, est/ref)

ccmat <- as.matrix(SetIndex(cc, "meas"))
bx <- barplot(ccmat, beside=TRUE, border=NA, axisnames=FALSE, ylim=c(0, 1.05), legend=TRUE)
xpos <- apply(bx, 2, mean)
xrg <- apply(bx, 2, range)
axis(1, bx, FALSE, tck=0.01)
text(xpos, par("usr")[3]-par("cxy")[2]*.2, adj=c(1,.5), colnames(ccmat), srt=45, xpd=NA)
box()
abline(h=c(1/1:4, 0), lty=2)

ggplot(nCfit)+
  geom_line(aes(time, ratio, color=meas))+
  lims(y=1+.2*c(-1,1))

write.csv(cc, "outputs/lambdaC_values_tseries.csv", row.names=FALSE)
