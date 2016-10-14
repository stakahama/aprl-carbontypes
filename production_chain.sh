#!/bin/sh

## -----------------------------------------------------------------------------

Rscript --vanilla select_SVOCs.R # SVOCs.json
Rscript --vanilla select_example_MCM.R # fig (apinene_ctype_tseries.pdf), fig (apinene_compound_abundance.pdf)

## -----------------------------------------------------------------------------

## clean directory
rm outputs/lambdaC_*

Rscript --vanilla production_lambdaC_compounds.R
Rscript --vanilla production_lambdaC_tseries.R
Rscript --vanilla production_lambdaC_coeff.R
Rscript --vanilla export_lambdaC.R
Rscript --vanilla production_ensemble.R

## -----------------------------------------------------------------------------

Rscript --vanilla figures_nC_grid.R # fig (fig_nCgrid_fg.pdf)
Rscript --vanilla figures_nC_cumsum.R # fig (nC_apinene_cumsum.pdf)
Rscript --vanilla figures_OSC_eval.R # fig (OSC_distr_meas.pdf)
Rscript --vanilla figures_props.R # fig (production_fig_props.pdf)
Rscript --vanilla figures_nC_reconstruction.R # fig (nC_est_scatterplot.pdf)
Rscript --vanilla figures_nC_tseries.R # fig (nC_recovery_tseries.pdf)
Rscript --vanilla figures_Sax.R # fig (fig_Sax.pdf)
Rscript --vanilla figures_nC_unc.R # fig (nC_unc_deltas.pdf)
