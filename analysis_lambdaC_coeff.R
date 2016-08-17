
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

df <- ldply(ff, function(f, x=1.96) {
  out <- readRDS(f)
  cc <- ResetIndex(as.data.frame(coef(summary(out))), "FG")
  cc[,"FG"] <- substring(cc[,"FG"], 2)
  cc[,c("bl", "bu")] <- cc[,"Estimate"] + outer(cc[,"Std. Error"], x*c(-1, 1), "*")
  cc
}, .id="meas")

levels(df$meas) <- with(list(x=levels(df$meas)),
                        c(sapply(split(x, ifelse(grepl("_collapsed", x), 1, 0)), sort)))

## -----------------------------------------------------------------------------

dodge <- position_dodge(width=0.9)
ggp <- ggplot(df)+
  geom_hline(yintercept=1/(1:3), linetype=2)+
  geom_bar(aes(FG, Estimate, fill=meas), stat="identity", position=dodge)+
  geom_errorbar(aes(FG, ymin=bu, ymax=bl, group=meas), width=0.025, position=dodge)+
  theme(axis.text.x=element_text(angle=30, hjust=1))

pdf(FilePath("lambdaC_coef_errorbars"), width=10, height=5)
print(ggp)
dev.off()

## -----------------------------------------------------------------------------0

nom <- 1/(1:3)
df$nominal <- nom[sapply(df$Estimate, function(x, y) which.min(abs(x-y)), nom)]

## with(df, plot(Estimate, nominal))
## abline(0, 1)

lambdaCarray <- acast(df, meas~FG, value.var="Estimate")
write.csv(lambdaCarray, FilePath("lambdaC_coef_actual"))

lambdaCarray <- acast(df, meas~FG, value.var="nominal")
write.csv(lambdaCarray, FilePath("lambdaC_coef_nominal"))
