#!/bin/sh


prefix=merged
ssp=~/git/projects/aprl-ssp
here=`pwd`

mkdir rev_data; cd $_

R -e "write.csv(unique(plyr::ldply(file.path('https://raw.githubusercontent.com/stakahama/aprl-ssp/master/validation',  c('apinenepropenemech.csv', 'tmbpropenemech.csv')), read.csv)), file='${prefix}_SMILES.csv', row.names=FALSE)"

# exit 1

python ${ssp}/substructure_search.py \
       -d -g common_atoms.csv \
       -i ${prefix}_SMILES.csv \
       -o ${prefix}_commonatoms.csv

python ${ssp}/substructure_search.py \
       -d -g SIMPOLgroups.csv -e SIMPOLexportlist.csv \
       -i ${prefix}_SMILES.csv \
       -o ${prefix}_SIMPOLgroups.csv

python ${ssp}/substructure_generate_fulltable.py \
       -g ${here}/patterns/MCMgroups.csv \
       -i ${prefix}_SMILES.csv \
       -o ${prefix}_MCMgroups

python ${ssp}/substructure_generate_fulltable.py \
       -d -g OScBonds.csv \
       -i ${prefix}_SMILES.csv \
       -o ${prefix}_OSc

python ${ssp}/substructure_adjacent_atoms.py \
       -i ${prefix}_SMILES.csv \
       -o ${prefix}_adjacent_atoms.csv

Rscript --vanilla ${ssp}/validation_atoms.R \
    ${prefix}_MCMgroups_atomfulltable.csv \
    ${prefix}_commonatoms.csv \
    ${prefix}

cd $here

