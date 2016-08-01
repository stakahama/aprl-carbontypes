library(limSolve)

## -----------------------------------------------------------------------------

## y is 1-D

lsq <- function(x, arg.A, arg.b)
  crossprod(arg.A %*% x - arg.b)

boundedls <- function(y, X, bl=NULL, bu=NULL, tol=.Machine$double.eps) {
  theta <- limSolve::nnls(X, y)
  init <- pmax(0L, pmin(theta$X, 1L))
  if(is.null(bl) && is.null(bu)) {
    DBind[G, H] <- buildconstr(X, bl=init-tol)
  } else {
    DBind[G, H] <- buildconstr(X, bl, bu)
  }
  out <- try(constrOptim(init, lsq, NULL, ui=G, ci=H, arg.A=X, arg.b=y), TRUE)
  if(inherits(out, "try-error")) browser()
  setNames(out$par, colnames(X))
}

buildconstr <- function(X, bl=0, bu=1, tol=with(.Machine, c(double.neg.eps, double.eps))) {
  if(length(bl)==1) bl <- rep(bl-tol[1], ncol(X))
  if(length(bu)==1) bu <- rep(bu+tol[2], ncol(X))
  G <- rbind(diag(ncol(X)), -diag(ncol(X)))
  H <- c(bl, -bu)
  list(G=G, H=H)
}

## -----------------------------------------------------------------------------

## Y is 2-D

lsq2 <- function(x, arg.A, arg.B, arg.dim) {
  norm(arg.B - arg.A %*% matrix(x, arg.dim[1], arg.dim[2]), "F")
}

buildconstr2 <- function(dim., bl=0, bu=1, tol=with(.Machine, c(double.neg.eps, double.eps))) {
  nrow. <- dim.[1]
  ncol. <- dim.[2]
  if (length(bl) == 1)
    bl <- rep(bl - tol[1], nrow.)
  if (length(bu) == 1)
    bu <- rep(bu + tol[2], nrow.)
  ##
  unitrow <- function(i, X) {
    X[i, ] <- 1
    X
  }
  zeros <- matrix(0, nrow., ncol.)
  ##
  V <- t(sapply(seq(nrow.), unitrow, zeros))
  G <- rbind(V, -V)
  H <- c(bl, -bu)
  ##
  list(G=G, H=H)
}

boundedls2 <- function (X, Y, bl = 0, bu = 1, init=0, ...) {
  dim. <- c(ncol(X), ncol(Y))
  dimnames. <- list(colnames(X), colnames(Y))
  n <- prod(dim.)
  ##
  ## *violates constraints
  ## seed <- 1
  ## set.seed(seed)
  ## sgn <- sample(c(-1, 1), n, TRUE)
  ## init <- sgn*rep(1/n, n)
  ##
  ## *violates constraints
  ## init <- solve(t(X) %*% X, t(X) %*% Y)
  ##
  if(identical(init, "n")) { ## produces more positive values than 0
    init <- rep(1/n, n)
  } else if(length(init)==1) {
    init <- rep(init, n)
  }
  ##
  DBind[G, H] <- buildconstr2(dim., bl, bu)
  out <- try(constrOptim(init, lsq2, NULL,
                         ui = G, ci = H,
                         ...,
                         arg.A = X,
                         arg.B = Y,
                         arg.dim = dim.),
             TRUE)
  if (inherits(out, "try-error"))
    browser()
  ##
  mat <- matrix(out[["par"]], dim.[1], dim.[2], dimnames=dimnames.)
  remaining <- setdiff(names(out), "par")
  attributes(mat)[remaining] <- out[remaining]
  ##
  mat
}

## test buildconstr2()
##
## m <- cbind(1:3, 1:3)
## dim. <- dim(m)
## r <- buildconstr2(dim.)
## all.equal(r$G %*% c(m), matrix(c(rowSums(m), -rowSums(m))))
