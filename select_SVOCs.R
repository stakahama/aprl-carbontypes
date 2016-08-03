
options(stringsAsFactors=FALSE)

library(dplyr)
library(Rfunctools)
library(RJSONIO)
PopulateEnv("mylib", "lib/lib_collapse.R")

## -----------------------------------------------------------------------------

AddPrefix <- function(x, prefix="merged")
  paste(prefix, x, sep="_")

inpfiles <- c(
  "molecattr"=file.path("data", AddPrefix("molec_attributes.csv")),
  "matrices"=file.path("data", AddPrefix("matrices.rda"))
)

outfiles <- c(
  "svoc"=file.path("inputs", "SVOCs.json"),
  "clabels"=file.path("inputs", "clabels.json")
)

## -----------------------------------------------------------------------------

DBind[X, Y, Theta, gamma] <- ReadRDA(inpfiles["matrices"])

molec.attr <- read.csv(inpfiles["molecattr"])

## -----------------------------------------------------------------------------

svoc <- setdiff(with(arrange(molec.attr, logC0), compound[round(logC0) <= 3]),
                Radicalgroups(X))

k <- apply(Y[svoc,] > 0, 2, any)
j <- apply(Theta[k,] > 0, 2, any)

## -----------------------------------------------------------------------------

out <- list(compounds=svoc, ctypes=rownames(Theta)[k], groups=colnames(Theta)[j])

clabels <- with(list(sk=out$ctypes, rest=setdiff(rownames(Theta), out$ctypes)),
                setNames(c(sprintf("S%02d",seq_along(sk)),
                           sprintf("T%02d",seq_along(rest))),
                         c(sk, rest)))

## Y <- Y[svoc,sk]
## X <- X[svoc,sj]
## Theta <- Theta[k,sj]
## gamma <- gamma[sj]

## -----------------------------------------------------------------------------

cat(toJSON(out), file=outfiles["svoc"])

cat(toJSON(clabels), file=outfiles["clabels"])
