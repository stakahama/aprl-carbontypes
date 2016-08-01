library(plyr)
library(dplyr)

Shorttype <- function(x)
  toupper(substring(x, 1, 1))

AtomOSc <- function(fulltable, ...) {
  ## OSc fulltable

  ## From Ruggeri and Takahama 2016
  z <- c(
    "C-H"=-1,
    "C-C"=0,
    "C=C"=0,
    "C-N"=1,
    "C-O"=1,
    "C=O"=2
  )

  ## Select carbon atoms and assign Ox numbers
  fulltable <- fulltable %>%
    filter(Shorttype(type)=="C") %>%
    mutate(z=z[as.character(group)])

  ## Carbon oxidation states
  fulltable <- fulltable %>%
    group_by(compound, atom) %>%
    summarize(OSc=sum(z)) %>%
    ungroup()

  ## Return
  fulltable

}

MolecOSc <- function(fulltable, ...) {

  fulltable <- AtomOSc(fulltable)

  ## Compound mean carbon oxidation states
  fulltable <- fulltable %>%
    group_by(compound) %>%
    summarize(OSc=mean(OSc))

  ## Return
  fulltable

}
