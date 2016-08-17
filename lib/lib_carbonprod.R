
CarbonInnerProd <- function(M, Y) {
  cmpds <- intersect(names(M), rownames(Y))
  time <- index(M)
  moles <- M[,cmpds] %*% Y[cmpds,]
  ## frac <- moles/rowSums(moles)
  melt(data.frame(time, moles, check.names=FALSE),
       id.vars="time",
       variable.name="clabel",
       value.name="nC")
}

CarbonProd <- function(M, Y) {
  cmpds <- intersect(names(M), rownames(Y))
  timevec <- index(M)
  M <- coredata(M)
  rownames(M) <- names(timevec) <- sprintf("x%d", seq_along(timevec))
  M.lf <- melt(M[,cmpds,drop=FALSE], varnames=c("time", "compound"), value.name="nmolec")
  Y.lf <- melt(Y[cmpds,], varnames=c("compound", "clabel"), value.name="nC") %>%
    mutate(clabel=factor(clabel, sort(levels(clabel))))
  full_join(M.lf, Y.lf, by=c("compound")) %>%
    group_by(time) %>% mutate(nC=nC*nmolec, nmolec=NULL) %>%
    ungroup() %>% mutate(time=timevec[time])
}

OrderSlice <- function(df) {
  ## add order to compound and carbon types
  wf <- acast(df, compound~clabel, sum, value.var="nC")
  levs <- rownames(wf)[order(rowSums(wf), decreasing=TRUE)]
  df %>% mutate(compound=factor(compound, levs),
                clabel=factor(clabel, sort(levels(clabel))))
}