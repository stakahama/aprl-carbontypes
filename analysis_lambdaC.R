
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("mylib", c("lib/lib_collapse.R", "lib/lib_constrOptim.R"))

## -----------------------------------------------------------------------------

AddPrefix <- function(x, prefix="merged")
  paste(prefix, x, sep="_")

OutFile <- function(x, ext=NULL, path="outputs", prefix="lambdaC")
  function(suffix=NULL)
    file.path(path, paste(sep=".", paste(sep="_", prefix, x, suffix), ext))

inpfiles <- c(
  "molecattr"=file.path("data", AddPrefix("molec_attributes.csv")),
  "matrices"=file.path("data", AddPrefix("matrices.rda")),
  "meas"=file.path("inputs", "meas_FGs.json"),
  "svoc"=file.path("inputs", "SVOCs.json")
)

outfiles <- list(
  "plot_ctype"=OutFile("ctype", "pdf"),
  "plot_nC"=OutFile("nC", "pdf"),
  "Phi"=OutFile("Phi", "rds"),
  "values"=OutFile("values", "csv"),
  "lmfit"=OutFile("lmfit", "rds")
)

## -----------------------------------------------------------------------------

CountMethod <- function(ni=1, nC, Theta, Y, gamma) {
  if(length(ni)==1)
    ni <- setNames(rep(ni, nrow(Y)), rownames(Y))
  1/(gamma*sum(ni*nC)) * (ni %*% Y %*% (1/rowSums(Theta)))
}

## -----------------------------------------------------------------------------

matrices <- ReadRDA(inpfiles["matrices"])

molec.attr <- read.csv(inpfiles["molecattr"])

measlist <- fromJSON(inpfiles["meas"])$lambdaC
collapserule <- fromJSON(inpfiles["meas"])$collapse

svoc <- fromJSON(inpfiles["svoc"])$compounds

## -----------------------------------------------------------------------------

DBind[measlist.collapsed, matrices.collapsed] <-
  AggGroups(collapserule, measlist, matrices)

## -----------------------------------------------------------------------------

loopvars <- list(c("matrices", "measlist"),
                 c("matrices.collapsed", "measlist.collapsed"))

for(elem in loopvars) {

  DBind[X, Y, Theta, gamma] <- get(elem[1])
  measlist <- get(elem[2])

  for(.label in names(measlist)) {

    meas <- measlist[[.label]]

    measC <- rowSums(Theta[,meas]) > 0
    Theta.s <- sweep(Theta[measC, meas], 2, gamma[meas],`*`)
    nC.s <- rowSums(sweep(Y[svoc,measC], 2, sign(rowSums(Theta.s)), `*`))
    Y.s <- Y[svoc,measC]
    X.s <- X[svoc,meas]
    gamma.s <- gamma[meas]

    ## -----------------------------------------------------------------------------

    ## direct calculation of Phi

    Phi <- list()

    Phi$solve <- solve(t(X.s) %*% X.s, t(X.s) %*% Y.s)

    Phi$ginv <- MASS::ginv(Theta.s)
    dimnames(Phi$ginv) <- rev(dimnames(Theta.s))

    Phi$constr <- BoundedLsq2(X.s, Y.s, init=0, outer.iterations=1e3)

    ## -----------------------------------------------------------------------------

    ## calculation of lambda

    lambdaC <- list()
    lambdaC$count <- CountMethod(1, nC.s, Theta.s, Y.s, gamma.s)
    lambdaC$solve <- solve(t(X.s) %*% X.s, t(X.s) %*% rowSums(Y.s))[,1]
    lambdaC$constr <- BoundedLsq(rowSums(Y.s), X.s, 0, 1)
    ## lambdaC$solve.rowsum <- rowSums(Phi$solve)
    lambdaC$constr.rowsum <- rowSums(Phi$constr)
    lambdaC$ginv.rowsum <- rowSums(Phi$ginv)

    ## -----------------------------------------------------------------------------

    ## C-type evaluation

    lf.ctype <- full_join(
      melt(Y.s, varnames=c("compound", "ctype"), value.name="ref"),
      ldply(Phi, function(Phi, X)
        melt(X %*% Phi, varnames=c("compound", "ctype"), value.name="est"),
        X.s, .id="method")
    )

    limits <- with(lf.ctype, data.frame(lim=range(c(ref, est))))

    ggp <- ggplot(lf.ctype)+
      geom_point(aes(ref, est, color=ctype), position=position_jitter(w=.1, h=0))+
      facet_wrap(~method)+
      geom_abline(intercept=0, slope=1)+
      guides(color=FALSE)+
      geom_blank(aes(lim, lim), data=limits)

    pdf(outfiles[["plot_ctype"]](.label), width=12, height=5)
    print(ggp)
    dev.off()

    ## -----------------------------------------------------------------------------

    ## nC evaluation

    lambdaC.mat <- do.call(rbind, lambdaC)

    lf.nC <- melt(data.frame(compound=rownames(Y.s), ref=rowSums(Y.s), X.s %*% t(lambdaC.mat)),
                  id.vars=c("compound", "ref"), variable.name="method", value.name="est")
    lf.nC <- left_join(lf.nC, select(molec.attr, compound, logC0))

    limits <- with(lf.nC, data.frame(lim=range(c(ref, est))))

    ggp <- ggplot(lf.nC)+
      geom_point(aes(ref, est, color=logC0), position=position_jitter(w=.25, h=0))+
      facet_wrap(~method)+
      geom_abline(intercept=0, slope=1)+
      geom_blank(aes(lim, lim), data=limits)

    pdf(outfiles[["plot_nC"]](.label))
    print(ggp)
    dev.off()

    ## -----------------------------------------------------------------------------

    ## export values

    saveRDS(Phi, file=outfiles[["Phi"]](.label))

    write.csv(ResetIndex(as.data.frame(lambdaC.mat), "method"),
              outfiles[["values"]](.label), row.names=FALSE)


    ## -----------------------------------------------------------------------------

    ## lm for confidence intervals

    fit.lm <- lm(y~X-1, data=list(y=rowSums(Y.s), X=X.s))
    saveRDS(fit.lm, file=outfiles[["lmfit"]](.label))

  }

}
