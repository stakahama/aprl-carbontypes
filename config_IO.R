
FilePath <- function(f) {
  filelist <- c(
    ##
    "mcmmass"="data-raw/^mcm_.+_mass\\.txt$",
    "simpol"="data/merged_SIMPOLgroups.csv",
    "fulltable"="data/merged_MCMGroups_atomfulltable.csv",
    "adjacent"="data/merged_adjacent_atoms.csv",
    "carbonattr"="data/merged_C_attributes.csv",
    "groupattr"="data/merged_group_attributes.csv",
    "molecattr"="data/merged_molec_attributes.csv",
    "matrices"="data/merged_matrices.rda",
    "matrices_2"="data/merged_matrices_2.rda",
    ##
    "svoc"="inputs/SVOCs.json",
    "clabels"="inputs/clabels.json",
    "meas"="inputs/meas_FGs.json",
    "example_1"="inputs/example_1.json",
    ## lambdaC
    "plot_ctype_fit"="outputs/lambdaC_ctype_%s.pdf",
    "plot_nC_fit"="outputs/lambdaC_nC_%s.pdf",
    "Phi"="outputs/lambdaC_Phi_%s.rds",
    "lambdaC"="outputs/lambdaC_values_%s.csv",
    "lmfit"="outputs/lambdaC_lmfit_%s.rds",
    ## OSc
    "plot_OSc"="outputs/OSc_OSc_%s_%s.pdf",
    "plot_elemratios"="outputs/OSc_elemratios_%s_%s.pdf",
    "plot_vankrevelen"="outputs/OSc_vanKrevelen_%s_%s.pdf",
    ## examples
    "tseries_gas"="inputs/apinene_formatted.csv",
    "tseries_aer"="inputs/apinene_aer_formatted.csv",
    "plot_ctype_tseries"="outputs/apinene_ctype_tseries.pdf",
    "plot_OSc_tseries"="outputs/apinene_OSc_tseries.pdf",
    "plot_compound_abundance"="outputs/apinene_compound_abundance.pdf",
    ##
    "plot_compound_nC"="outputs/nC_compound.pdf",
    "plot_nC_cumsum"="outputs/nC_apinene_cumsum.pdf"
  )
  filelist[f]
}

GenericReader <- function(filename) {
  ext <- toupper(sub(".+\\.(.+)$", "\\1", basename(filename)))
  Read <- switch(
    ext,
    "JSON"=RJSONIO::fromJSON,
    "CSV"=read.csv,
    "RDS"=readRDS,
    "RDA"=Rfunctools::ReadRDA
  )
  Read(filename)
}

ReadFile <- function(f, ...) {
  Read <- switch(f,
                 "simpol"=function(x) as.matrix(read.csv(x, row.names=1)),
                 GenericReader)
  Read(FilePath(f))
}


SprintF <- function(f, ...)
  sprintf(FilePath(f), ...)

## OutFile <- function(x, ext=NULL, path="outputs", prefix="lambdaC")
##   function(suffix=NULL)
##     file.path(path, paste(sep=".", paste(sep="_", prefix, x, suffix), ext))
