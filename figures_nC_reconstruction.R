options(stringsAsFactors=FALSE)

library(abind)
library(plyr)
library(dplyr)
library(reshape2)
library(Rfunctools)
library(RJSONIO)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_collapse.R", "lib/lib_lambdaC.R"))
GGTheme()

## -----------------------------------------------------------------------------

matrices <- ReadFile("matrices")
DBind[measlist, collapserule] <-
  ReadFile("meas")[c("lambdaC", "collapse")]
svoc <- ReadFile("svoc")$compounds
molec.attr <- ReadFile("molecattr")
decisions <- as.list(ReadFile("example_1"))

lambdaC <- ReadFile("lambdaC_array")

## -----------------------------------------------------------------------------
DBind[X, Y, Theta, gamma] <- matrices

methods <- c("count", "solve", "fit", "nominal")

out <- list()
for(.meas in names(measlist)) {
  j <- intersect(measlist[[.meas]], colnames(Theta))
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
levels(out$meas) <- Capitalize(levels(out$meas))

table <- left_join(out,
                   molec.attr %>%
                   mutate(logC0=round(logC0)) %>%
                   select(compound, logC0, OSC))

## table <- out

Stats <- function(x) {
  out <- lm(est~ref-1, data=x)
  val <- c(slope=unname(coef(out)), cor=with(x, cor(ref, est)))
  data.frame(as.list(val))
}

## table %>% count(meas, method)
stats <- table %>% group_by(meas, method) %>%
  do(Stats(.))

stats.formatted <- stats
stats.formatted$slope <- sprintf("ratio = %.2f", stats$slope)
stats.formatted$cor <- sprintf("italic(r) == %.2f", stats$cor)

grid <- do.call(expand.grid, lapply(out[c("method", "meas")], levels))
grid$letter <- sprintf("%s)", letters[seq(nrow(grid))])

yoffset <- 1.25
ggp <- ggplot(table)+
  geom_point(aes(ref, est, color=OSC))+
  facet_grid(meas~method)+
  geom_abline(intercept=0, slope=1)+
  scale_x_continuous(name=expression(italic(n)[C]^"*"), limits=c(2, 11), breaks=seq(2, 11, 2))+
  scale_y_continuous(name=expression(hat(italic(n))[C]^"*"), limits=c(2, 11), breaks=seq(2, 11, 2))+
  geom_text(aes(x=-Inf, y=Inf, label=letter), data=grid, size=5, hjust=0, vjust=1)+
  ## scale_color_continuous(name=expression(logC[0]))+
  scale_color_continuous(name=expression(bar(OS)[C]), low=colors.OSC[1], high=tail(colors.OSC, 1))+
  geom_text(aes(x=-Inf, y=Inf, label=slope), data=stats.formatted, size=5, hjust=0, vjust=1+yoffset)+
  geom_text(aes(x=-Inf, y=Inf, label=cor), data=stats.formatted, size=5, hjust=0, vjust=1+2*yoffset, parse=TRUE)

pdf(FilePath("plot_nC_est"), width=9, height=6)
print(ggp)
dev.off()
