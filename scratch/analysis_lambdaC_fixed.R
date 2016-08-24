
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", c("lib/lib_collapse.R", "lib/lib_constrOptim.R", "lib/lib_lambdaC.R"))

## -----------------------------------------------------------------------------

matrices <- ReadFile("matrices")
molec.attr <- ReadFile("molecattr")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
svoc <- ReadFile("svoc")$compounds

## -----------------------------------------------------------------------------

uniq <- with(matrices, UniqueMapping(Theta))

DBind[measlist.collapsed, matrices.collapsed, uniq.collapsed] <-
  AggGroups(collapserule, measlist, matrices, uniq)

## -----------------------------------------------------------------------------

matrices.original <- matrices
measlist.original <- measlist
uniq.original <- uniq

loopvars <- list(c("matrices.original", "measlist.original", "uniq.original"),
                 c("matrices.collapsed", "measlist.collapsed", "uniq.collapsed"))

for(.elem in loopvars) {

  DBind[X, Y, Theta, gamma] <- get(.elem[1])
  measlist <- get(.elem[2])
  uniq <- get(.elem[3])

  multi <- with(subset(uniq, group %in% group[duplicated(group)]),
                split(data.frame(ctype, value), group))

  if(length(multi) > 0) {
    for(i in names(multi)) {
      oldlab <- multi[[i]]$ctype
      mval <- multi[[i]]$value
      newlab <- paste(oldlab, collapse="+")
      ##
      Y <- cbind(Y[,setdiff(colnames(Y),oldlab)],
                 last=Y[,oldlab] %*% mval)
      colnames(Y)[ncol(Y)] <- newlab
      ##
      Theta <- rbind(Theta[setdiff(rownames(Theta),oldlab),],
                     Theta[oldlab[1],])
      rownames(Theta)[nrow(Theta)] <- newlab
      ##
      uniq <- rbind(subset(uniq, !ctype %in% oldlab),
                    data.frame(ctype=newlab, group=i, value=mval[1]))
    }
  }

  for(.label in names(measlist)) {

    meas <- measlist[[.label]]
    measC <- names(which(rowSums(Theta[,meas]) > 0))

    fullm <- list(ctype=measC, group=meas)
    extra <- list(ctype=intersect(uniq$ctype, measC),
                  group=intersect(uniq$group, meas))

    Y.u <- Y[svoc,extra$ctype,drop=FALSE]
    X.u <- X[svoc,extra$group,drop=FALSE]

    meas <- setdiff(meas, uniq$group)
    measC <- setdiff(measC, uniq$ctype)
    Theta.s <- sweep(Theta[measC,meas], 2, gamma[meas],`*`)
    nC.s <- rowSums(Y[svoc,measC])
    Y.s <- Y[svoc,measC]
    X.s <- X[svoc,meas]
    gamma.s <- gamma[meas]

    ## -----------------------------------------------------------------------------

    ExpandPhi <- function(Phi) {
      n <- sapply(fullm[c("group", "ctype")], length)
      p <- matrix(0, n[1], n[2], dimnames=fullm[c("group", "ctype")])
      grid <- do.call(expand.grid, c(dimnames(Phi), list(stringsAsFactors=FALSE)))
      for(i in seq(nrow(grid)))
        p[grid[i,1],grid[i,2]] <- Phi[grid[i,1],grid[i,2]]
      for(i in seq(nrow(uniq)))
        if(uniq[i,"group"] %in% rownames(p) && uniq[i,"ctype"] %in% colnames(p))
          p[uniq[i,"group"],uniq[i,"ctype"]] <- uniq[i,"value"]
      p
    }

    ExpandLambdaC <- function(lambdaC) {
      lam <- with(fullm, setNames(numeric(length(group)), group))
      lam[names(lambdaC)] <- lambdaC
      for(i in seq(nrow(uniq)))
        if(uniq[i,"group"] %in% names(lam))
          lam[uniq[i,"group"]] <- uniq[i,"value"]
      lam
    }

    ## -----------------------------------------------------------------------------

    ## direct calculation of Phi

    Phi <- list()

    Phi$solve <- ApproxSolve(X.s, Y.s)

    Phi$ginv <- MASS::ginv(Theta.s)
    dimnames(Phi$ginv) <- rev(dimnames(Theta.s))

    Phi$constr <- BoundedLsq2(X.s, Y.s, init=0, outer.iterations=1e3)

    Phi[] <- lapply(Phi, ExpandPhi)

    ## -----------------------------------------------------------------------------

    ## calculation of lambda

    lambdaC <- list()
    ## lambdaC$count <- CountMethod(1, nC.s, Theta, Y[svoc,], gamma, meas)
    lambdaC$count <- CountMethod2(Theta, gamma, meas)
    lambdaC$solve <- ApproxSolve(X.s, Vec2Mat(rowSums(Y.s)))[,1]
    lambdaC$constr <- BoundedLsq(rowSums(Y.s), X.s, 0, 1)
    ## lambdaC$solve.rowsum <- rowSums(Phi$solve)
    lambdaC$constr.rowsum <- rowSums(Phi$constr)
    lambdaC$ginv.rowsum <- rowSums(Phi$ginv)

    lambdaC[] <- lapply(lambdaC, ExpandLambdaC)

    ## -----------------------------------------------------------------------------

    ## C-type evaluation

    lf.ctype <- full_join(
      melt(cbind(Y.s, Y.u), varnames=c("compound", "ctype"), value.name="ref"),
      ldply(Phi, function(Phi, X)
        melt(X[,rownames(Phi)] %*% Phi, varnames=c("compound", "ctype"), value.name="est"),
        cbind(X.s, X.u), .id="method")
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

    lf.nC <- melt(data.frame(compound=rownames(Y.s),
                             ref=rowSums(cbind(Y.s, Y.u)),
                             cbind(X.s, X.u)[,colnames(lambdaC.mat)] %*% t(lambdaC.mat)),
                  id.vars=c("compound", "ref"), variable.name="method", value.name="est")
    lf.nC <- left_join(lf.nC, select(molec.attr, compound, logC0))

    limits <- with(lf.nC, data.frame(lim=range(c(ref, est))))

    ggp <- ggplot(lf.nC)+
      geom_point(aes(ref, est, color=logC0), position=position_jitter(w=.1, h=0))+
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
