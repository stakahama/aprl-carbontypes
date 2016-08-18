
am <- c("C"=12.01, "H"=1.001, "N"=14.35, "O"=15.99)
heteroatoms <- c("H", "N", "O")

Calculate <- function(n, X, Y, Theta, gamma, zFG, Lambda, i=TRUE, k=TRUE, j=TRUE, nC.total=NULL) {

  ## extract
  Theta.s <- sweep(Theta[k,j], 2, gamma[j],`*`)
  nC.s <- n[i]*rowSums(sweep(Y[i,k], 2, sign(rowSums(Theta.s)), `*`))
  Y.s <- n[i]*Y[i,k]
  X.s <- n[i]*X[i,j]
  gamma.s <- gamma[j]
  zFG.s <- zFG[j]
  Lambda.s <- Lambda[c("H","N","O"), j]

  ## local variables
  ## scoped variable
  ## `am`

  ## *** molecule basis ***

  OC <- am["C"]*nC.s
  OSC <- (X.s %*% zFG.s)/nC.s
  atoms <- X.s %*% t(Lambda.s)
  ratios <- atoms/nC.s
  masses <- sweep(atoms, 2, am[colnames(atoms)], "*")
  OM <- OC+rowSums(masses)
  OM.OC <- OM/OC
  OM.OC.comp <- cbind(C=1, masses/OC)

  basis.molec <- NamedList(nC.total, ratios, OSC, OM, OM.OC, OM.OC.comp)

  ## *** mixture basis ***

  ## compute
  if(!is.numeric(nC.total))
    nC.total <- sum(nC.s)

  ##
  OC <- am["C"]*nC.total
  OSC <- sum(X.s %*% zFG.s)/nC.total
  atoms <- X.s %*% t(Lambda.s)
  ratios <- colSums(atoms)/nC.total
  masses <- am[colnames(atoms)]*colSums(atoms)
  OM <- unname(OC+sum(masses))
  OM.OC <- unname(OM/OC)
  OM.OC.comp <- c(C=1, masses/OC)

  basis.mix <- NamedList(nC.total, ratios, OSC, OM, OM.OC, OM.OC.comp)

  ## return
  list(molec=basis.molec, mix=basis.mix)

}

Calculate2 <- function(n, X, Y, Theta, gamma, zFG, Lambda, zeta, cOM, i=TRUE, k=TRUE, j=TRUE, lambdaC=NULL) {

  if(isTRUE(i)) i <- rownames(Y)
  if(isTRUE(k)) k <- rownames(Theta)
  if(isTRUE(j)) j <- colnames(Theta)

  ## extract
  Theta.s <- sweep(Theta[k,j], 2, gamma[j],`*`)
  nC.s <- n[i]*rowSums(sweep(Y[i,k], 2, sign(rowSums(Theta.s)), `*`))
  Y.s <- n[i]*Y[i,k]
  X.s <- n[i]*X[i,j]
  gamma.s <- gamma[j]
  zFG.s <- zFG[j]
  Lambda.s <- Lambda[heteroatoms, j] # heteroatoms lexically scoped
  zeta.s <- zeta[k]

  ## compute
  if(is.numeric(lambdaC)) {
    lambdaC.s <- lambdaC[j]
    lambdaC.s[] <- replace(lambdaC.s, is.na(lambdaC.s), 0)
    nC.s <- sum(X.s %*% lambdaC.s)
  }

  nC.total <- sum(nC.s)
  nC.c <- colSums(n[i]*sweep(Y[i,k], 2, sign(rowSums(Theta.s)), `*`))
  OC.c <- am["C"]*nC.c
  OM.c <- nC.c*cOM[names(nC.c)]
  atoms.g <- sweep(Lambda.s, 2, colSums(X.s), "*")
  atomr.g <- atoms.g/nC.total
  OM.OC.g <- colSums(am[rownames(atomr.g)]*atoms.g)/sum(OC.c)

  OSC.true <- sum(sweep(Y.s, 2, zeta.s, "*"))/nC.total
  ## OSC.approx <- sum(X.s %*% zFG.s)/nC.total
  OSC.approx <- (t(colSums(X.s)) %*% zFG.s)/nC.total

  ctype.subset <- names(which(nC.c > 0))
  group.subset <- names(which(OM.OC.g > 0))

  ## *** generate tables ***

  ## masses
  masses <- inner_join(data.frame(ctype=names(OC.c), OC=OC.c),
                       data.frame(ctype=names(OM.c), OM=OM.c)) %>%
    filter(ctype %in% ctype.subset) #%>% melt(., "ctype")

  ## ratios
  rownames(atomr.g) <- paste0(rownames(atomr.g), "/C")
  ratios <- full_join(adply(atomr.g, 2, .id="group"), #melt(atomr.g, c("metric", "group"))
                      data.frame(group=names(OM.OC.g),
                                 "OM/OC-1"=OM.OC.g,
                                 check.names=FALSE)) %>%
    filter(group %in% group.subset)

  ## OSC
  oxstate <- with(list(x=c(true=OSC.true, approx=OSC.approx)),
                  data.frame(method=names(x), value=x))

  ## return
  list(masses=masses, ratios=ratios, OSC=oxstate)

}

CarbonTypeMass <- function(carbon.attr) {
  ##
  id <- "ctype" #c("type", "ctype")
  ##
  carbon.mass <- unique(carbon.attr[,c(id, heteroatoms)])
  carbon.mass[heteroatoms] <- sweep(carbon.mass[heteroatoms], 2, am[heteroatoms], "*")
  carbon.mass[["C"]] <- am["C"]
  carbon.mass[["OM"]] <- carbon.mass[["C"]]+rowSums(carbon.mass[heteroatoms])
  ##
  carbon.mass

}
