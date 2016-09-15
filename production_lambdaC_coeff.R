
options(stringsAsFactors=FALSE)

library(Rfunctools)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")
PopulateEnv("fig", "config_fig.R")
PopulateEnv("mylib", c("lib/lib_collapse.R"))
GGTheme()

## -----------------------------------------------------------------------------

ff <- list.files("outputs", "lmfit", full=TRUE)
names(ff) <- gsub("^lambdaC\\_lmfit\\_(.+)\\.rds$","\\1",basename(ff))

f2 <- list.files("outputs", "lambdaC_count", full=TRUE)
names(f2) <- gsub("^lambdaC\\_count\\_(.+)\\.csv$","\\1",basename(f2))

f3 <- ReadFile("lambdaC_tseries")

Theta <- ReadFile("matrices")[["Theta"]]

measlist <- ReadFile("meas")[["lambdaC"]]
decisions <- as.list(ReadFile("example_1"))

## -----------------------------------------------------------------------------

df <- ldply(ff, function(f, x=1.96) {
  out <- readRDS(f)
  cc <- ResetIndex(as.data.frame(coef(summary(out))), "group")
  cc[,"group"] <- substring(cc[,"group"], 2)
  cc[,c("bl", "bu")] <- cc[,"Estimate"] + outer(cc[,"Std. Error"], x*c(-1, 1), "*")
  cc
}, .id="meas")

df2 <- ldply(f2, function(f, x=1.96) {
  cc <- ResetIndex(t(as.matrix(read.csv(f, row.names="stat", check.names=FALSE))), "group")
  cc[,c("bl", "bu")] <- cc[,"mean"] + outer(cc[,"se"], x*c(-1, 1), "*")
  cc
}, .id="meas")

df3 <- melt(f3, id.vars="meas", variable.name="group", value.name="Estimate")

uniq <- UniqueMapping(Theta)

## -----------------------------------------------------------------------------

Sortfn <- function(x)
  c(sapply(split(x, ifelse(grepl("_collapsed", x), 1, 0)), sort))

dfm <- full_join(df %>% mutate(method="solve"),
                 df2 %>% rename(Estimate=mean, `Std. Error`=se) %>%
                 mutate(method="count")) %>%
  mutate(meas=factor(meas, Sortfn(unique(meas))))

dfm <- full_join(dfm, df3 %>% mutate(method="fit"))

dfm <- dfm %>% filter(!group %in% uniq$group)

dodge <- position_dodge(width=0.9)
## ggp <- ggplot(dfm)+
##   facet_grid(method~.)+
##   geom_hline(yintercept=c(1/(1:4), 0), linetype=2)+
##   geom_bar(aes(group, Estimate, fill=meas), stat="identity", position=dodge)+
##   geom_errorbar(aes(group, ymin=bu, ymax=bl, group=meas), width=0.025, position=dodge)+
##   theme(axis.text.x=element_text(angle=30, hjust=1))

ggp <- local({
  methods <- c("count", "solve", "fit")
  dfm$method <- factor(dfm$method, methods, toupper(methods))
  dfm$meas <- with(dfm, factor(meas, unique(meas), Capitalize(unique(meas))))

  ggplot(dfm)+
    facet_grid(meas~.)+
    geom_hline(yintercept=c(1/(1:4), 0), linetype=2)+
    geom_bar(aes(group, Estimate, fill=method), stat="identity", position=dodge)+
    geom_errorbar(aes(group, ymin=bu, ymax=bl, group=method), width=0.025, position=dodge)+
    theme(axis.text.x=element_text(angle=30, hjust=1))+
    scale_y_continuous(limits=c(0, 1.8), expand=c(0, 0))+
    labs(x="", y=expression("Coefficient,"~lambda[C,j]))
})

pdf(FilePath("lambdaC_coef_errorbars"), width=7, height=7)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------


## fixed.values <- with(uniq, setNames(value, group))[fixed]

## lambdaCtable <- dcast(dfm, method+meas~group, value.var="Estimate")
## lambdaCtable <- cbind(lambdaCtable, fixed.values)
## lambdaCtable <- lambdaCtable[,c("method", "meas", intersect(colnames(Theta), names(lambdaCtable)))]
## write.csv(lambdaCtable, FilePath("lambdaC_coef_actual"), row.names=FALSE)

fixed.groups <- intersect(Reduce(union, as.list(unname(measlist))), uniq$group)
fixed.df <- with(uniq, data.frame(group, Estimate=value, fixed=TRUE)) %>%
  filter(group %in% fixed.groups)

dfm2 <- dfm %>% group_by(meas, method) %>%
  do(full_join(., cbind(.[1,c("meas", "method")], fixed.df)))

dfm2$group <- factor(dfm2$group, intersect(colnames(Theta), unique(dfm2$group)))
dfm2$fixed <- ifelse(is.na(dfm2$fixed), FALSE, dfm2$fixed)

saveRDS(dfm2, FilePath("lambdaC_coef_actual_f"))
