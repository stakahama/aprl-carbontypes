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

OutFile <- function(x, ext=NULL, path="outputs", prefix="OSc")
  function(suffix=NULL)
    file.path(path, paste(sep=".", paste(sep="_", prefix, x, suffix), ext))

inpfiles <- c(
  "matrices"=file.path("data", AddPrefix("matrices.rda")),
  "matrices_2"=file.path("data", AddPrefix("matrices_2.rda")),
  "molecattr"=file.path("data", AddPrefix("molec_attributes.csv")),
  "meas"=file.path("inputs", "meas_FGs.json"),
  "svoc"=file.path("inputs", "SVOCs.json")
)

outfiles <- list(
  "plot_compound_OSc"=OutFile("OSc", "pdf"),
  "plot_compound_elemratios"=OutFile("elemratios", "pdf"),
  "plot_compound_vankrevelen"=OutFile("vanKrevelen", "pdf")
  ## "Phi"=OutFile("Phi", "rds"),
  ## "values"=OutFile("values", "csv"),
  ## "lmfit"=OutFile("lmfit", "rds")
)

## -----------------------------------------------------------------------------

matrices <- c(ReadRDA(inpfiles["matrices"]), ReadRDA(inpfiles["matrices_2"]))

molec.attr <- read.csv(inpfiles["molecattr"])

measlist <- fromJSON(inpfiles["meas"])$lambdaC
collapserule <- fromJSON(inpfiles["meas"])$collapse

svoc <- fromJSON(inpfiles["svoc"])$compounds

## -----------------------------------------------------------------------------

DBind[measlist.collapsed, matrices.collapsed] <-
  AggGroups(collapserule, measlist, matrices[c("X","Y","Theta","gamma","zFG", "Lambda")])

## -----------------------------------------------------------------------------

nC <- rowSums(matrices$Y[svoc,])
atoms <- with(matrices, X[svoc,] %*% t(Lambda[c("H","N","O"),]))

## -----------------------------------------------------------------------------

loopvars <- list(c("matrices", "measlist"),
                 c("matrices.collapsed", "measlist.collapsed"))

for(elem in loopvars) {

  DBind[X, Y, Theta, gamma, zFG, Lambda] <- get(elem[1])
  measlist <- get(elem[2])

  ## nC <- rowSums(Y)

  for(.label in names(measlist)) {

    meas <- measlist[[.label]]

    measC <- rowSums(Theta[,meas]) > 0
    Theta.s <- sweep(Theta[measC, meas], 2, gamma[meas],`*`)
    nC.s <- rowSums(sweep(Y[svoc,measC], 2, sign(rowSums(Theta.s)), `*`))
    Y.s <- Y[svoc,measC]
    X.s <- X[svoc,meas]
    gamma.s <- gamma[meas]
    zFG.s <- zFG[meas]

    ## -------------------------------------------------------------------------

    ## *** compound OSc ***

    paired <- data.frame(compound=svoc,
                         OSc=(X[svoc,] %*% zFG)/nC[svoc],
                         OSc.s=(X.s %*% zFG.s)/nC.s)
    paired <- left_join(paired, select(molec.attr, compound, logC0))

    ggp <- ggplot(paired)+
      geom_point(aes(OSc, OSc.s, color=logC0))+
      lims(x=c(-3, 2), y=c(-3, 2))+
      labs(x=expression(OS[C]), y=expression(OS[C]*"*"))+
      geom_abline(intercept=0, slope=1)

    pdf(outfiles[["plot_compound_OSc"]](.label))
    print(ggp)
    dev.off()

    ## -------------------------------------------------------------------------

    ## *** compound element ratios ***

    ## atoms <- X[svoc,] %*% t(Lambda[c("H","N","O"),])
    atoms.s <- X.s %*% t(Lambda[c("H","N","O"),meas])

    paired <- full_join(
      melt(atoms, varnames=c("compound", "elem"), value.name="full"),
      melt(atoms.s, varnames=c("compound", "elem"), value.name="star")
    )

    limits <- paired %>% group_by(elem) %>%
      do(with(., data.frame(value=c(full, star))))

    ggp <- ggplot(paired)+
      geom_point(aes(full, star))+
      geom_blank(aes(value, value), limits)+
      facet_wrap(~elem, scale="free")+
      geom_abline(intercept=0, slope=1)

    pdf(outfiles[["plot_compound_elemratios"]](.label), width=12, height=5)
    print(ggp)
    dev.off()

    ## -------------------------------------------------------------------------

    ## *** compound van Krevelen diagrams ***

    lf <- melt(paired %>%
               mutate(full=full/nC[as.character(compound)],
                      star=star/nC.s[as.character(compound)]),
               measure.vars=c("full", "star"))

    wf <- dcast(lf, compound+variable~elem)
    wwf <- dcast(lf, compound~elem+variable)

    ggp <- ggplot(wf)+
      geom_segment(aes(H_full, O_full, xend=H_star, yend=O_star),
                   arrow=arrow(length = unit(0.015, "npc")),
                   data=wwf, color=8)+
      geom_point(aes(H, O, color=variable), shape=0)+
      lims(x=c(1, 3), y=c(0, 2))+
      labs(x=expression(H/C), y=expression(O/C))

    pdf(outfiles[["plot_compound_vankrevelen"]](.label))
    print(ggp)
    dev.off()

    ## break()

  }

}
