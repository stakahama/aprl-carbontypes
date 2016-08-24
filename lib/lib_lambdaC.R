
CountMethod <- function(ni=1, nC, Theta, Y, gamma, J.s) {
  if(length(ni)==1)
    ni <- setNames(rep(ni, nrow(Y)), rownames(Y))
  ## 1/(gamma*sum(ni*nC)) * (ni %*% Y %*% (1/rowSums(Theta)))
  nC.s <- nC
  lambda.C <- setNames(numeric(length(J.s)), J.s)
  for(j in J.s) {
    k <- Theta[,j] > 0
    ## gamma.s <- gamma[j]
    Theta.s <- Theta[k,,drop=FALSE]
    Y.s <- Y[,k,drop=FALSE]
    lambda.C[j] <- 1/(sum(ni %*% Y.s)) *
      (ni %*% Y.s %*% (1/(Theta.s %*% gamma)))
  }
  lambda.C
}

## CountMethod2 <- function(Theta, gamma, J.s) {
##   lambda.C <- setNames(numeric(length(J.s)), J.s)
##   for(j in J.s) {
##     k <- Theta[,j] > 0
##     Theta.s <- Theta[k,,drop=FALSE]
##     lambda.C[j] <- mean(1/(Theta.s %*% gamma))
##   }
##   lambda.C
## }

CountMethod2 <- function(Theta, gamma, J.s) {
  lambdaC <- setNames(numeric(length(J.s)), J.s)
  n <- lambdaC.se <- lambdaC
  for(j in J.s) {
    k <- Theta[,j] > 0
    Theta.s <- Theta[k,,drop=FALSE]
    reciprocal <- 1/(Theta.s %*% gamma)
    lambdaC[j] <- mean(reciprocal)
    n[j] <- length(which(k))
    lambdaC.se[j] <- sd(reciprocal)/sqrt(n[j])
  }
  list(lambdaC, lambdaC.se, n)
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


Nominal <- function(lambdaC, method="mean", adjust=NULL) {
  nom <- c(1/sort(unique(rowSums(Theta))), 0)
  template <- apply(lambdaC, c(1,3), mean)
  if(method=="mean") {
    mat <- template
  } else {
    mat <- lambdaC[,method,]
  }
  template[] <- sapply(mat, function(x, y) y[which.min(abs(x-y))], nom)
  if(!is.null(adjust))
    for(j in names(adjust))
      template[,j] <- adjust[j]
  template
}
