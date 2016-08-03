
Radicalgroups <- function(X)
  rownames(X)[apply(X[,grepl("radical", colnames(X))] > 0, 1, any)]

AggGroups <- function(agg, meas, matrices) {

  newmatrices <- matrices
  newmeas <- list()
  for(new in names(agg)) {
    ##
    DBind[X, Y, Theta, gamma, zFG, Lambda] <- newmatrices
    ##
    target <- agg[[new]]
    groups <- colnames(Theta)
    rest <- setdiff(groups, target)
    ##
    Theta <- cbind(Theta[,rest,drop=FALSE], rowSums(Theta[,target,drop=FALSE]))
    colnames(Theta) <- c(rest, new)
    ##
    gamma <- c(gamma[rest], mean(gamma[target]))
    names(gamma) <- c(rest, new)
    ##
    if(!is.null(zFG)) {
      ## this approximation may be problematic
      zFG <- c(zFG[rest], mean(zFG[target]))
      names(zFG) <- c(rest, new)
    }
    if(!is.null(Lambda)) {
      ## this approximation may be problematic
      Lambda <- cbind(Lambda[,rest,drop=FALSE], rowMeans(Lambda[,target,drop=FALSE]))
      colnames(Lambda) <- c(rest, new)
    }
    ##
    X <- sweep(Y %*% Theta, 2, gamma, `*`)
    ##
    newmatrices <- list(X=X, Y=Y, Theta=Theta, gamma=gamma, zFG=zFG, Lambda=Lambda)
    ##
  }
  for(set in names(meas)) {
    newset <- paste(sep="_", set, "collapsed")
    newvars <- meas[[set]]
    for(new in names(agg))
      newvars <- c(setdiff(newvars, agg[[new]]), new)
    newmeas[[newset]] <- intersect(colnames(X), newvars)
  }

  list(meas=newmeas, matrices=newmatrices)

}
