
FilePath <- function(f) {
  filelist <- c(
    ## -----------------------------------------------------------------------------
    ## "mcmmass"="data-raw/^mcm_.+_mass\\.txt$",
    ## "simpol"="data/merged_SIMPOLgroups.csv",
    ## "fulltable"="data/merged_MCMGroups_atomfulltable.csv",
    ## "adjacent"="data/merged_adjacent_atoms.csv",
    ## "carbonattr"="data/merged_C_attributes.csv",
    ## "groupattr"="data/merged_group_attributes.csv",
    ## "molecattr"="data/merged_molec_attributes.csv",
    ## "matrices"="data/merged_matrices.rda",
    ## "matrices_2"="data/merged_matrices_2.rda",
    ## -----------------------------------------------------------------------------
    "mcmmass"="data-raw/^mcm_apinene_mass\\.txt$",
    "simpol"="data/apinene_SIMPOLgroups.csv",
    "fulltable"="data/apinene_MCMGroups_atomfulltable.csv",
    "adjacent"="data/apinene_adjacent_atoms.csv",
    "carbonattr"="data/apinene_C_attributes.csv",
    "groupattr"="data/apinene_group_attributes.csv",
    "molecattr"="data/apinene_molec_attributes.csv",
    "matrices"="data/apinene_matrices.rda",
    "matrices_2"="data/apinene_matrices_2.rda",
    ## -----------------------------------------------------------------------------
    "svoc"="inputs/SVOCs.json",
    "clabels"="inputs/clabels.json",
    "meas"="inputs/meas_FGs.json",
    "example_1"="inputs/example_1.json",
    ## -----------------------------------------------------------------------------
    ## lambdaC
    "plot_ctype_fit"="outputs/lambdaC_ctype_%s.pdf",
    "plot_nC_fit"="outputs/lambdaC_nC_%s.pdf",
    "Phi"="outputs/lambdaC_Phi_%s.rds",
    "lambdaC"="outputs/lambdaC_values_%s.csv",
    "lambdaC_count"="outputs/lambdaC_count_%s.csv",
    "lmfit"="outputs/lambdaC_lmfit_%s.rds",
    ## OSC
    "plot_OSC"="outputs/OSC_OSC_%s_%s.pdf",
    "plot_elemratios"="outputs/OSC_elemratios_%s_%s.pdf",
    "plot_vankrevelen"="outputs/OSC_vanKrevelen_%s_%s.pdf",
    ## examples
    "tseries_gas"="inputs/apinene_formatted.csv",
    "tseries_aer"="inputs/apinene_aer_formatted.csv",
    "plot_ctype_tseries"="outputs/apinene_ctype_tseries.pdf",
    "plot_OSC_tseries"="outputs/apinene_OSC_tseries.pdf",
    "plot_compound_abundance"="outputs/apinene_compound_abundance.pdf",
    ## nC
    "plot_compound_nC"="outputs/nC_compound.pdf",
    "plot_nC_cumsum"="outputs/nC_apinene_cumsum.pdf",
    ## nC (2)
    "plot_nC_tseries"="outputs/nC_tseries_estimated.pdf",
    ## properties
    "plot_props"="outputs/props_%s.pdf",
    "plot_props_scatter"="outputs/props_scatter_%s.pdf",
    "props_file"="outputs/propsfile_%s.rds",
    "plot_props_cumsum"="outputs/props_cumsum_%s.pdf",
    "tables_props_cumsum"="outputs/props_cumsum_%s.rds",
    ## -----------------------------------------------------------------------------
    ## lambdaC (f)
    "lambdaC_coef_errorbars"="outputs/lambdaC_coef_errorbars.pdf",
    "lambdaC_coef_nominal"="outputs/lambdaC_coef_nominal.csv",
    "lambdaC_coef_actual"="outputs/lambdaC_coef_actual.csv",
    "lambdaC_coef_actual_f"="outputs/lambdaC_coef_actual.rds",
    "lambdaC_tseries"="outputs/lambdaC_values_tseries.csv",
    "lambdaC_array"="outputs/lambdaC_array.rds",
    ## properties (f)
    "plot_props_1"="outputs/production_fig_props_1.pdf",
    "plot_props_2"="outputs/production_fig_props_2.pdf",
    "plot_osc_f"="outputs/production_fig_OSC.pdf",
    "plot_osc_distr"="outputs/OSC_distr_meas.pdf",
    ## nC (f)
    "plot_nC_est"="outputs/nC_est_scatterplot.pdf",
    "plot_nC_recovery"="outputs/nC_recovery_tseries.pdf",
    ## nC unc (f)
    "plot_nC_unc"="outputs/nC_unc_deltas.pdf",
    ## Sax (f)
    "plot_Sax"="outputs/fig_Sax.pdf",
    ## nC_grid (f)
    "plot_nCgrid_ctype"="outputs/fig_nCgrid_ctype.pdf",
    "plot_nCgrid_fg"="outputs/fig_nCgrid_fg.pdf",
    ##
    "lastline"=NA
  )
  filelist[f]
}

GenericReader <- function(filename) {
  ext <- toupper(sub(".+\\.(.+)$", "\\1", basename(filename)))
  Read <- switch(
    ext,
    "JSON"=RJSONIO::fromJSON,
    "CSV"=function(...) read.csv(..., check.names=FALSE),
    "RDS"=readRDS,
    "RDA"=Rfunctools::ReadRDA
  )
  Read(filename)
}

ReadFile <- function(f, ...) {
  if(file.exists(f)) {
    GenericReader(f)
  } else {
    Read <- switch(f,
                   "simpol"=MatrixReader,
                   ## "lambdaC_coef_nominal"=MatrixReader,
                   ## "lambdaC_coef_actual"=MatrixReader,
                   GenericReader)
    Read(FilePath(f))
  }
}

MatrixReader <- function(x)
  as.matrix(read.csv(x, row.names=1, check.names=FALSE))

SprintF <- function(f, ...)
  sprintf(FilePath(f), ...)

## OutFile <- function(x, ext=NULL, path="outputs", prefix="lambdaC")
##   function(suffix=NULL)
##     file.path(path, paste(sep=".", paste(sep="_", prefix, x, suffix), ext))
