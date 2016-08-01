library(plyr)
library(dplyr)

Shorttype <- function(x)
  toupper(substring(x, 1, 1))

Groupvars <- function(fulltable)
  sort(unique(fulltable$group))

WidenGroupsC <- function(fulltable)
  reshape2::dcast(fulltable %>% filter(Shorttype(type) == "C" & is.finite(match)),
                  compound + atom + type ~ group, length, value.var = "match")

AddCtypeWide <- function(fulltable) {
  wf <- WidenGroupsC(fulltable)
  fgvars <- Groupvars(fulltable)
  wf$ctype <- apply(wf[, fgvars], 1, Int2Key)
  wf
}

Int2Key <- function(x)
  MultiIndex(x, collapse=",", format=pryr::partial(sprintf, "%d"))
