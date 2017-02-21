#!/bin/sh

# --- do manually ---
# cp -p config_IO.R config_IO_apinene.R
# cp -p rev_config_IO.R config_IO.R

Rscript --vanilla build_molec_attributes.R
Rscript --vanilla build_group_attributes.R
Rscript --vanilla build_matrices.R
Rscript --vanilla build_matrices_2.R
Rscript --vanilla build_C_attributes.R

# --- do manually ---
# mv config_IO_apinene.R config_IO.R
