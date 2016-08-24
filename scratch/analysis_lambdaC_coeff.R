
options(stringsAsFactors=FALSE)

library(Rfunctools)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())
PopulateEnv("IO", "config_IO.R")

## -----------------------------------------------------------------------------

ff <- list.files("outputs", "lmfit", full=TRUE)
names(ff) <- gsub("^lambdaC\\_lmfit\\_(.+)\\.rds$","\\1",basename(ff))

f2 <- list.files("outputs", "lambdaC_count", full=TRUE)
names(f2) <- gsub("^lambdaC\\_count\\_(.+)\\.csv$","\\1",basename(f2))

Theta <- ReadFile("matrices")[["Theta"]]

## -----------------------------------------------------------------------------

df <- ldply(ff, function(f, x=1.96) {
  out <- readRDS(f)
  cc <- ResetIndex(as.data.frame(coef(summary(out))), "FG")
  cc[,"FG"] <- substring(cc[,"FG"], 2)
  cc[,c("bl", "bu")] <- cc[,"Estimate"] + outer(cc[,"Std. Error"], x*c(-1, 1), "*")
  cc
}, .id="meas")

df2 <- ldply(f2, function(f, x=1.96) {
  cc <- ResetIndex(t(as.matrix(read.csv(f, row.names="stat", check.names=FALSE))), "FG")
  cc[,c("bl", "bu")] <- cc[,"mean"] + outer(cc[,"se"], x*c(-1, 1), "*")
  cc
}, .id="meas")


## -----------------------------------------------------------------------------

Sortfn <- function(x)
  c(sapply(split(x, ifelse(grepl("_collapsed", x), 1, 0)), sort))

dfm <- full_join(df %>% mutate(method="solve"),
                 df2 %>% rename(Estimate=mean, `Std. Error`=se) %>%
                 mutate(method="count")) %>%
  mutate(meas=factor(meas, Sortfn(unique(meas))))

dodge <- position_dodge(width=0.9)
ggp <- ggplot(dfm)+
  facet_grid(method~.)+
  geom_hline(yintercept=c(1/(1:4), 0), linetype=2)+
  geom_bar(aes(FG, Estimate, fill=meas), stat="identity", position=dodge)+
  geom_errorbar(aes(FG, ymin=bu, ymax=bl, group=meas), width=0.025, position=dodge)+
  theme(axis.text.x=element_text(angle=30, hjust=1))

pdf(FilePath("lambdaC_coef_errorbars"), width=10, height=5)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------0

## dfm.agg <- dfm %>% group_by(meas, FG) %>% summarize(Estimate=mean(Estimate))

nom <- c(1/sort(unique(rowSums(Theta))), 0)
dfm$nominal <- nom[sapply(dfm$Estimate, function(x, y) which.min(abs(x-y)), nom)]

## with(df, plot(Estimate, nominal))
## abline(0, 1)

## lambdaCarray <- acast(dfm.agg, meas~FG, value.var="Estimate")
## write.csv(lambdaCarray, FilePath("lambdaC_coef_actual"))

## lambdaCarray <- acast(dfm.agg, meas~FG, value.var="nominal")
## write.csv(lambdaCarray, FilePath("lambdaC_coef_nominal"))


lambdaCtable <- dcast(dfm, method+meas~FG, value.var="Estimate")
write.csv(lambdaCtable, FilePath("lambdaC_coef_actual"), row.names=FALSE)

lambdaCtable <- dcast(dfm, method+meas~FG, value.var="nominal")
write.csv(lambdaCtable, FilePath("lambdaC_coef_nominal"), row.names=FALSE)
