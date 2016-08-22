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
PopulateEnv("mylib", "lib/lib_collapse.R")
GGTheme()

## -----------------------------------------------------------------------------

matrices <- ReadFile("matrices")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
svoc <- ReadFile("svoc")$compounds
molec.attr <- ReadFile("molecattr")

## patt <- "lambdaC_Phi_(.+_collapsed).rds"
## ff <- list.files("outputs", patt, full=TRUE)
## names(ff) <- sub(patt, "\\1", basename(ff))

ff <- list.files("outputs", "lambdaC_values.+_collapsed\\.csv$", full=TRUE)
names(ff) <- gsub("(lambdaC_values_)|(\\.csv)", "", basename(ff))

## -----------------------------------------------------------------------------

DBind[measlist.collapsed, matrices.collapsed] <-
  AggGroups(collapserule, measlist, matrices)

DBind[X, Y, Theta, gamma] <- matrices.collapsed

## -----------------------------------------------------------------------------


lambdaC <- ldply(ff, read.csv, check.names=FALSE, .id="meas") %>%
  filter(method %in% c("count", "solve"))

ReplZero <- function(x) replace(x, is.na(x), 0)
id.vars <- c("meas", "method")

lambdaC <- daply(lambdaC, id.vars, function(x, i)
  ReplZero(unlist(x[setdiff(names(x), i)])), id.vars)

methods <- c("count", "solve")

out <- list()
for(.meas in names(measlist.collapsed)) {
  j <- measlist.collapsed[[.meas]]
  nC.ref <- rowSums(sweep(Y[svoc,], 2, sign(Theta[,j] %*% gamma[j]), "*"))
  nC.est <- list()
  for(.method in methods) {
    value <- X[svoc,j] %*% lambdaC[.meas,.method,j]
    nC.est[[.method]] <- data.frame(compound=svoc, ref=nC.ref, est=value)
  }
  nC.est <- ldply(nC.est, .id="method")
  out[[.meas]] <- nC.est
}
out <- ldply(out, .id="meas")

levels(out$method) <- toupper(levels(out$method))
levels(out$meas) <- Capitalize(sub("_collapsed", "", levels(out$meas)))

table <- left_join(out,
                   molec.attr %>%
                   mutate(logC0=round(logC0)) %>%
                   select(compound, logC0, OSC))

## table <- out

grid <- do.call(expand.grid, lapply(out[c("method", "meas")], levels))
grid$letter <- sprintf("%s)", letters[seq(nrow(grid))])

ggp <- ggplot(table)+
  geom_point(aes(ref, est, color=OSC))+
  facet_grid(meas~method)+
  geom_abline(intercept=0, slope=1)+
  scale_x_continuous(name=expression(italic(n)[C]^"*"), limits=c(2, 11), breaks=seq(2, 11, 2))+
  scale_y_continuous(name=expression(hat(italic(n))[C]^"*"), limits=c(2, 11), breaks=seq(2, 11, 2))+
  geom_text(aes(x=-Inf, y=Inf, label=letter), data=grid, size=5, hjust=0, vjust=1)+
  ## scale_color_continuous(name=expression(logC[0]))+
  scale_color_continuous(name=expression(bar(OS)[C]), low=colors.OSC[1], high=tail(colors.OSC, 1))

pdf("outputs/nC_est_scatterplot.pdf", width=5, height=7)
print(ggp)
dev.off()
