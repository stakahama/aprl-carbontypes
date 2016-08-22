# APRL Carbontypes

Code to accompany manuscript (*in prep*) on carbon type analysis.

## Repository structure

Folders:

* lib/: code library for reused functions
* data-raw/: original data
* data/: populated with build\_ssp.sh
* outputs/: output of analysis\_*

Scripts:

* build\_*: data prep (for constructing master tables). read/write to data/
* analysis\_*: analysis. write to outputs/
* select\_*: selection of example
* config\_: configurations
* figures\_*: production figures


## Scripts

1. build\_ssp.sh
2. build\_attributes.sh


## Notes

- This is the working directory.
- `options(stringsAsFactors=FALSE)` is the default.
- "config_IO.R" facilitates I/O across files

