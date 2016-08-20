
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_collapse.R", "lib/lib_constrOptim.R"))

## -----------------------------------------------------------------------------

CountMethod <- function(ni=1, nC, Theta, Y, gamma) {
  if(length(ni)==1)
    ni <- setNames(rep(ni, nrow(Y)), rownames(Y))
  1/(gamma*sum(ni*nC)) * (ni %*% Y %*% (1/rowSums(Theta)))
}

## (need for phenol in set4)
ApproxSolve <- function(X, Y, B=matrix(0, ncol(X), ncol(Y), dimnames=list(colnames(X), colnames(Y))), thresh=1, inc=1) {
  ## reduce dimensionality of the problem but keep the dimensions
  ##   of the original solution (if it worked) in the output
  out <- try(solve(t(X) %*% X, t(X) %*% Y), TRUE)
  if(class(out)=="try-error") {
    j <- names(which(colSums(X) > thresh))
    ApproxSolve(X[,j], Y, B, thresh+inc) ## recursive!
  } else {
    B[rownames(out), colnames(out)] <- out
    B
  }
}

Vec2Mat <- function(x, column="Y") {
  x <- as.matrix(x)
  colnames(x) <- column
  x
}

## -----------------------------------------------------------------------------

matrices <- ReadFile("matrices")
molec.attr <- ReadFile("molecattr")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
svoc <- ReadFile("svoc")$compounds

## -----------------------------------------------------------------------------

DBind[measlist.collapsed, matrices.collapsed] <-
  AggGroups(collapserule, measlist, matrices)

## -----------------------------------------------------------------------------

matrices.original <- matrices
measlist.original <- measlist

loopvars <- list(c("matrices.original", "measlist.original"),
                 c("matrices.collapsed", "measlist.collapsed"))

for(.elem in loopvars) {

  DBind[X, Y, Theta, gamma] <- get(.elem[1])
  measlist <- get(.elem[2])

  for(.label in names(measlist)) {

    meas <- measlist[[.label]]

    ## thresh <- 1 # > 0 produces singularity
    ## meas <- intersect(meas, names(which(colSums(X[svoc,meas]) > thresh)))

    measC <- rowSums(Theta[,meas]) > 0
    Theta.s <- sweep(Theta[measC, meas], 2, gamma[meas],`*`)
    nC.s <- rowSums(sweep(Y[svoc,measC], 2, sign(rowSums(Theta.s)), `*`))
    Y.s <- Y[svoc,measC]
    X.s <- X[svoc,meas]
    gamma.s <- gamma[meas]

    ## -----------------------------------------------------------------------------

    ## direct calculation of Phi

    Phi <- list()

    Phi$solve <- ApproxSolve(X.s, Y.s)

    Phi$ginv <- MASS::ginv(Theta.s)
    dimnames(Phi$ginv) <- rev(dimnames(Theta.s))

    Phi$constr <- BoundedLsq2(X.s, Y.s, init=0, outer.iterations=1e3)

    ## -----------------------------------------------------------------------------

    ## calculation of lambda

    lambdaC <- list()
    lambdaC$count <- CountMethod(1, nC.s, Theta.s, Y.s, gamma.s)
    lambdaC$solve <- ApproxSolve(X.s, Vec2Mat(rowSums(Y.s)))[,1]
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

    pdf(SprintF("plot_ctype_fit", .label), width=12, height=5)
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

    pdf(SprintF("plot_nC_fit", .label))
    print(ggp)
    dev.off()

    ## -----------------------------------------------------------------------------

    ## export values

    saveRDS(Phi, file=SprintF("Phi", .label))

    write.csv(ResetIndex(as.data.frame(lambdaC.mat), "method"),
              SprintF("lambdaC", .label), row.names=FALSE)

    ## -----------------------------------------------------------------------------

    ## lm for confidence intervals

    fit.lm <- lm(y~X-1, data=list(y=rowSums(Y.s), X=X.s))
    saveRDS(fit.lm, file=SprintF("lmfit", .label))

  }

}
