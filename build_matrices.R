
################################################################################
##
## build_matrices.R
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Oct. 2016
##
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################


options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)
PopulateEnv("IO", "config_IO.R")
PopulateEnv("mylib", "lib/lib_C_attributes.R")

## -----------------------------------------------------------------------------

fulltable <- ReadFile("fulltable")

## -----------------------------------------------------------------------------

Uniquify <- compose(unique, subset)

## -----------------------------------------------------------------------------

wf.groups <- AddCtypeWide(fulltable)
fgvars <- Groupvars(fulltable)

df.Y <- dcast(wf.groups %>% distinct(compound, atom, ctype),
              compound ~ ctype,
              length, value.var="atom")

df.Theta <- Uniquify(wf.groups,,c("ctype", fgvars))

## -----------------------------------------------------------------------------

df.X <- dcast(Uniquify(fulltable,, c(compound, match, group)),
              compound ~ group,
              length, value.var="match")

df.gamma <- Uniquify(fulltable %>% filter(Shorttype(type)=="C") %>%
                     count(compound, match, group) %>% ungroup(),,
                     c(group, n))

## -----------------------------------------------------------------------------

X <- as.matrix(SetIndex(df.X, "compound"))
Y <- as.matrix(SetIndex(df.Y, "compound"))
Theta <- as.matrix(SetIndex(df.Theta, "ctype"))[colnames(Y),,drop=FALSE]
gamma <- with(df.gamma, setNames(1/n, group))[colnames(X)]

## -----------------------------------------------------------------------------

save(X, Y, Theta, gamma, file=FilePath("matrices"))
