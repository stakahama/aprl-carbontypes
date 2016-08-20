options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_collapse.R", "lib/lib_constrOptim.R"))

## -----------------------------------------------------------------------------

load(FilePath("matrices"))
carbon.attr <- ReadFile("carbonattr")

## -----------------------------------------------------------------------------

LabelC <- function(x) {

  f <- function(count) {
    group <- names(count)
    sprintf("%s%s", ifelse(count==1, "", sprintf("%d",count)), group)
  }

  count <- x[x > 0]
  lab <- f(count)
  if(length(lab) > 1) { ## polygroup
    lab <- MultiIndex(lab)
  }

  lab

}

## colnames(Theta) <- Relabel(colnames(Theta), labels.FG)
## apply(Theta, 1, LabelC)

## -----------------------------------------------------------------------------

id <- c("ctype","type")
atomlist <- c("C","H","N","O")
groups <- colnames(X)

ctable <- unique(subset(carbon.attr,,c(id,atomlist,groups)))
ctable$index <- with(ctable, paste(type, ctype))

desig.C <- c("none", "primary", "secondary", "tertiary", "quaternary")
ctable$subst <- with(ctable, ifelse(type=="C3", desig.C[as.integer(C)+1L], NA))
ctable$lab <- apply(setNames(ctable[,groups], Relabel(groups, labels.FG)), 1, LabelC)

ctable[,c("type", "subst", "lab")]

clabels <- with(unique(ctable[,c("ctype", "lab")]), setNames(lab, ctype))

cat(toJSON(clabels), file=FilePath("clabels"))


## ctable %>% filter(type=="C3" & alcohol == 1 & `alkane CH`==2 & C == 1)
## ctable %>% filter(type=="C3" & alcohol == 1 & `alkane CH`==1 & C == 2)
## ctable %>% filter(type=="C3" & alcohol == 1 & `alkane CH`==0 & C == 3)

