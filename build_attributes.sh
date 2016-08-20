#!/bin/sh

Rscript --vanilla build_molec_attributes.R
Rscript --vanilla build_group_attributes.R
Rscript --vanilla build_matrices.R
Rscript --vanilla build_matrices_2.R
Rscript --vanilla build_C_attributes.R
