#!/bin/sh

## this is necessary
rm outputs/lambdaC_*

Rscript --vanilla select_SVOCs.R
Rscript --vanilla select_MCM.R # fig 1 (apinene_ctype_tseries.pdf), fig 2 (apinene_compound_abundance.pdf)

# Rscript --vanilla figures_nC_cumsum.R # fig 3 (nC_apinene_cumsum.pdf)

Rscript --vanilla production_lambdaC_compounds.R
Rscript --vanilla production_lambdaC_tseries.R
Rscript --vanilla production_lambdaC_coeff.R
Rscript --vanilla production_ensemble.R

Rscript --vanilla figures_props.R # fig (production_fig_props.pdf)
Rscript --vanilla figures_OSC_eval.R # fig (OSC_distr_meas.pdf)

Rscript --vanilla production_nC_tseries.R # fig (nC_recovery_tseries.pdf)
Rscript --vanilla figures_nC_reconstruction.R # fig (nC_est_scatterplot.pdf)
