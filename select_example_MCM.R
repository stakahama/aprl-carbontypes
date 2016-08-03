

options(stringsAsFactors=FALSE)

library(zoo)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(pryr)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("mylib", c("lib/lib_io.R", "lib/lib_collapse.R"))

## -----------------------------------------------------------------------------

AddPrefix <- function(x, prefix="merged")
  paste(prefix, x, sep="_")

Example <- partial(AddPrefix, prefix="apinene")

inpfiles <- c(
  "tseries_gas"=file.path("inputs", Example("formatted.csv")),
  "tseries_aer"=file.path("inputs", Example("aer_formatted.csv")),
  "carbonattr"=file.path("data", AddPrefix("C_attributes.csv")),
  "svoc"=file.path("inputs", "SVOCs.json"),
  "clabels"=file.path("inputs", "clabels.json"),
  "matrices"=file.path("data", AddPrefix("matrices.rda"))
)

outfiles <- c(
  "plot_ctype_tseries"=file.path("outputs", Example("ctype_tseries.pdf")),
  "plot_OSc_tseries"=file.path("outputs", Example("OSc_tseries.pdf")),
  "plot_compound_abundance"=file.path("outputs", Example("compound_abundance.pdf"))
)

## -----------------------------------------------------------------------------

DBind[svoc, sk, sj] <- fromJSON(inpfiles["svoc"])

DBind[X, Y, Theta, gamma] <- ReadRDA(inpfiles["matrices"])

moles.molec <- lapply(setNames(inpfiles[c("tseries_gas", "tseries_aer")], c("gas", "aer")),
                      ReadTSeries)

clabels <- fromJSON(inpfiles["clabels"])

carbon.attr <- read.csv(inpfiles["carbonattr"])

## -----------------------------------------------------------------------------

colnames(Y) <- clabels[rownames(Theta)]
rownames(Theta) <- clabels[rownames(Theta)]

carbon.attr$clabel <- clabels[carbon.attr$ctype]

## -----------------------------------------------------------------------------

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
  wf <- acast(df, compound~clabel, sum, value.var="nC")
  levs <- rownames(wf)[order(rowSums(wf), decreasing=TRUE)]
  df %>% mutate(compound=factor(compound, levs),
                clabel=factor(clabel, sort(levels(clabel))))
}

## -----------------------------------------------------------------------------

moles.lf <- ldply(moles.molec, CarbonInnerProd, Y, .id="phase")

frac.lf <- moles.lf %>% group_by(phase, time) %>%
  mutate(frac=nC/sum(nC))

ggp <- ggplot(frac.lf)+
  geom_area(aes(time, frac, fill=clabel))+
  facet_grid(phase~.)+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(expand=c(0, 0))

pdf(outfiles["plot_ctype_tseries"], width=10, height=7)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

osc.lf <- left_join(moles.lf, carbon.attr %>% select(clabel, OSc), by="clabel") %>%
  group_by(phase, time, OSc) %>% summarize(nC=sum(nC)) %>% ungroup() %>%
  group_by(phase, time) %>% mutate(frac=nC/sum(nC)) %>% ungroup()

## ggp <- ggplot(osc.lf %>% filter(time %in% seq(0, 25, 1)))+
##   geom_bar(aes(time, frac, fill=OSc), stat="identity")+
##   facet_grid(phase~.)+
##   scale_x_continuous(expand=c(0, 0))+
##   scale_y_continuous(expand=c(0, 0))
## print(ggp)

levs <- sort(unique(osc.lf$OSc), decreasing=TRUE)
colorscale <- colorRampPalette(rev(c("#132B43", "#56B1F7")))(length(levs))

ggp <- ggplot(osc.lf %>% mutate(OSc=factor(OSc, levs)))+
  geom_area(aes(time, frac, fill=OSc), color="white", size=0)+
  scale_fill_manual(values=colorscale)+
  facet_grid(phase~.)+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(expand=c(0, 0))

pdf(outfiles["plot_OSc_tseries"], width=10, height=7)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

decisions <- list(
  hour=20,
  ncompounds=20
)

example <- Slice(moles.molec$aer, decisions$hour)

example <- CarbonProd(example, Y)

example <- OrderSlice(example)

example.subset <- example %>%
  filter(unclass(compound) < decisions$ncompounds) %>%
  mutate(clabel=factor(clabel))

ggp <- ggplot(example.subset) +
  geom_bar(aes(compound, nC, fill=clabel), stat="identity")+
  theme(axis.text.x=element_text(angle=60, hjust=1))

pdf(outfiles["plot_compound_abundance"], width=10, height=7)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------

