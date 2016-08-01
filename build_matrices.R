
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(pryr)
PopulateEnv("mylib", "lib/lib_C_attributes.R")

File <- function(x, path="data", prefix="merged")
  file.path(path, paste(prefix, x, sep="_"))

## -----------------------------------------------------------------------------

inpfiles <- c(
  "fulltable"=File("MCMGroups_atomfulltable.csv")
)

outfiles <- c(
  "matrices"=File("matrices.rda")
)

## -----------------------------------------------------------------------------

fulltable <- read.csv(inpfiles["fulltable"])

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

save(X, Y, Theta, gamma, file=outfiles["matrices"])
