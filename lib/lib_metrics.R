
am <- c("C"=12.01, "H"=1.001, "N"=14.35, "O"=15.99)

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

Calculate2 <- function(n, X, Y, Theta, gamma, zFG, Lambda, carbon.attr, i=TRUE, k=TRUE, j=TRUE, lambdaC=NULL) {

  ## extract
  Theta.s <- sweep(Theta[k,j], 2, gamma[j],`*`)
  nC.s <- n[i]*rowSums(sweep(Y[i,k], 2, sign(rowSums(Theta.s)), `*`))
  Y.s <- n[i]*Y[i,k]
  X.s <- n[i]*X[i,j]
  gamma.s <- gamma[j]
  zFG.s <- zFG[j]
  Lambda.s <- Lambda[c("H","N","O"), j]
  lambdaC.s <- lambdaC[j]
  lambdaC.s[] <- replace(lambdaC.s, is.na(lambdaC.s), 0)

  carbon.mass <- CarbonTypeMass(carbon.attr)

  ## compute
  if(!is.numeric(lambdaC)) {
    nC.total <- sum(nC.s)
  } else {
    nC.total <- sum(X.s %*% lambdaC.s)
  }

  nC.c <- colSums(n[i]*sweep(Y[i,k], 2, sign(rowSums(Theta.s)), `*`))
  OC.c <- am["C"]*nC.c
  OM.c <- nC.c*with(carbon.mass, setNames(OM, ctype))[names(nC.c)]
  atoms.g <- sweep(Lambda.s, 2, colSums(X.s), "*")
  atomr.g <- atoms.g/nC.total
  OM.OC.g <- colSums(am[rownames(atomr.g)]*atoms.g)/sum(OC.c)

  ctype.subset <- names(which(nC.c > 0))
  group.subset <- names(which(OM.OC.g > 0))

  ## generate tables
  masses <- inner_join(data.frame(ctype=names(OC.c), OC=OC.c),
                       data.frame(ctype=names(OM.c), OM=OM.c)) %>%
    filter(ctype %in% ctype.subset) #%>% melt(., "ctype")

  rownames(atomr.g) <- paste0(rownames(atomr.g), "/C")
  ratios <- full_join(adply(atomr.g, 2, .id="group"), #melt(atomr.g, c("metric", "group"))
                      data.frame(group=names(OM.OC.g),
                                 "OM/OC-1"=OM.OC.g,
                                 check.names=FALSE)) %>%
    filter(group %in% group.subset)

  ## return
  list(masses=masses, ratios=ratios)

}

CarbonTypeMass <- function(carbon.attr) {
  ##
  id <- "ctype" #c("type", "ctype")
  atomlist <- c("C", "H", "N", "O")
  ##
  carbon.mass <- unique(carbon.attr[,c(id, atomlist)])
  carbon.mass$C <- 1
  carbon.mass[atomlist] <- sweep(carbon.mass[atomlist], 2, am[atomlist], "*")
  carbon.mass$OM <- rowSums(carbon.mass[atomlist])
  ##
  carbon.mass
}
