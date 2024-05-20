
## Lightcurve_Pipeline
A Pipeline written mostly by Erica Thygesen and Jessica Ranshaw designed to generate lightcurves from TESS and Kepler data that are flattened, chopped, and in a format readable by Jason Eastman's EXOFASTv2 code.

Requires Andrew Vanderburg's keplersplinev2 to run: https://github.com/avanderburg/keplersplinev2.

## get_FeH.ipynb

A script to generate a metallicity prior from the listed spectroscopic metallicities on ExoFOP or from a TRES data file.

## remake_pdf.ipynb

Remakes the posterior distribution functions (PDFs) from EXOFAST outputs for fits that are bimodal in stellar mass and age or for the radii of planets in grazing transit configurations.

## RV_mini_pipeline.ipynb

A very simple and incomplete script that takes a TRES or CHIRON data file and generates an RV file in the format EXOFAST requires for fits.

## get_dilution.ipynb

A short script that generates a dilution prior for a source observed by TESS by querying the TICv8.2 catalog using Vizier.

## generate_priors_TOI.ipynb

A script that generates a spectroscopic metallicity prior, TESS dilution prior, and planetary ephemeris starting points for a given TESS Object of Interest.
