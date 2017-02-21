
################################################################################
##
## rev_analysis.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################


options(stringsAsFactors=FALSE)

library(xtable)
library(Rfunctools)
library(plyr)
library(tidyverse)
PopulateEnv("lib", "lib/lib_units.R")
PopulateEnv("IO_apin", "config_IO.R")
PopulateEnv("IO", "rev_config_IO.R")
PopulateEnv("fig", "config_fig.R")

## -----------------------------------------------------------------------------

DBind[, Y, Theta, gamma, zFG, ] <-
  c(ReadFile("matrices"), ReadFile("matrices_2"))

colnames(Theta) <- Relabel(colnames(Theta), labels.FG)
names(gamma) <- Relabel(names(gamma), labels.FG)
names(zFG) <- Relabel(names(zFG), labels.FG)

onlyCgrp <- colnames(Theta) %in% c("tertiary sp2 carbon", "aromatic sp2 carbon", "quaternary carbon")
skelCgrp <- gamma < 1

counts <- sort(setNames(colSums(Y), colnames(Y)), decreasing=TRUE)
labels <- setNames(sprintf("X%d", seq_along(counts)), names(counts))

## -----------------------------------------------------------------------------

sum0p <- rowSums(Theta[,onlyCgrp]) > 0
sum1p <- rowSums(Theta[,!onlyCgrp & !skelCgrp]) > 0
sum2p <- rowSums(Theta[,!onlyCgrp & skelCgrp]) > 0

## classif <- ifelse(sum0p, "unfunctionalized",
##            ifelse(sum1p & !sum2p, "singlecarbon",
##            ifelse(sum1p & sum2p, "mixed",
##            ifelse(!sum1p & sum2p, "skeletalcarbon", NA))))

classif <- ifelse(sum0p, "unfunctionalized",
           ifelse(sum1p & !sum2p, "singlecarbon",
           ifelse(sum2p, "skeletalcarbon", NA)))

table(classif)

## -----------------------------------------------------------------------------

Boolvec <- function(x)
  classif %in% x

RemZero <- function(x) {
  theta <- Theta[Boolvec(x),]
  jj <- colSums(theta) > 0
  theta[,jj]
}

Mat2DF <- function(x) {
  zFG. <- zFG[colnames(x)]
  jj <- order(zFG., names(zFG.))
  OSC <- rowSums(sweep(x, 2, zFG., "*"))
  ##
  df <- data.frame(
    Label = labels[rownames(x)],
    OSC = OSC,
    n = counts[rownames(x)],
    `colnames<-`(x[,jj,drop=FALSE],
                 sprintf("\\rot{%s [%d]}", names(zFG.)[jj], zFG.[jj])),
    row.names=NULL,
    check.names=FALSE
  )
  df %>% arrange(OSC, desc(n))
}

Formatdf <- function(df) {
  ## continued from Mat2DF
  df$OSC <- sprintf("%.0f", df$OSC)
  df$n <- sprintf("%d", df$n)
  names(df)[2:3] <- c("OS$_\\chem{C}$", "$n$")
  ##
  df
}

## -----------------------------------------------------------------------------

subtables <- lapply(setNames(,unique(classif)), RemZero)

t(sapply(subtables, dim))

sum(Y)

for(x in names(subtables)) {
  sink(file.path("rev_outputs", sprintf("%s_table.tex", x)))
  df <- Mat2DF(subtables[[x]])
  xt <- xtable(Formatdf(df))
  align(xt) <- paste0("ll", paste(rep("c", ncol(df)-1), collapse=""))
  print(xt,
        include.rownames=FALSE,
        sanitize.colnames.function = identity)
  sink()
}

tbl <- ldply(subtables, Mat2DF, .id="type")

stopifnot(sum(tbl$n) == sum(Y))

tbl %>% group_by(type) %>% summarize(n=sum(n))

## -----------------------------------------------------------------------------

DBind[, Y.apin, Theta.apin, gamma.apin, zFG.apin, ] <-
  c(get("ReadFile", "IO_apin")("matrices"), get("ReadFile", "IO_apin")("matrices_2"))

colnames(Theta.apin) <- Relabel(colnames(Theta.apin), labels.FG)
names(gamma.apin) <- Relabel(names(gamma.apin), labels.FG)
names(zFG.apin) <- Relabel(names(zFG.apin), labels.FG)

newFG <- setdiff(colnames(Theta), colnames(Theta.apin))
carbonFG <- c("tertiary sp2 carbon", "quaternary carbon")

indexr <- apply(Theta[,colnames(Theta.apin)], 1, MultiIndex, collapse=",")
newp <- rowSums(Theta[,newFG]) > 0
carbonp <- rowSums(Theta[,carbonFG]) > 0
samep <- (indexr %in% rownames(Theta.apin)) & !newp
singleCp.apin <- samep & !carbonp

## length(which(samep & !newp))
## indexr[!samep & !newp]
## table(indexr)[setdiff(indexr, rownames(Theta.apin))]

singleCp <- !sum0p & sum1p & !sum2p

setdiff(rownames(Theta.apin), unique(indexr)) ## 0
rownames(Theta)[!samep] ## 19 new
rownames(Theta)[singleCp] ## 46 single C (combine)
rownames(Theta)[singleCp.apin] ## 39 single C (apinene)

singleC.orig <- rownames(Theta)[singleCp.apin]
singleC.new <- rownames(Theta)[singleCp & !singleCp.apin]

labx <- setNames(names(labels), labels)


## single carbon only
tbl %>%
  mutate(labx = labx[Label],
         class = ifelse(labx %in% singleC.orig, "orig",
                 ifelse(labx %in% singleC.new, "new",
                        "other"))) %>%
  group_by(class) %>%
  summarize(n = sum(n))

## new vs. old (41 types)
tbl %>%
  mutate(labx = labx[Label],
         class = ifelse(labx %in% rownames(Theta)[samep], "orig", "new")) %>%
  group_by(class) %>%
  summarize(n = sum(n))

## -----------------------------------------------------------------------------

clabels <- ReadFile("clabels")

df <- data.frame(old=clabels[indexr[samep]], new=labels[names(indexr)[samep]]) %>%
  arrange(old)

df$old <- sprintf("%d", df$old)

print(xtable(df), include.rownames=FALSE)

