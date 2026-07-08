
## Lightcurve_Pipeline
A Pipeline written mostly by Erica Thygesen and Jessica Ranshaw designed to generate lightcurves from TESS and Kepler data that are flattened, chopped, and in a format readable by Jason Eastman's EXOFASTv2 code.

Requires Andrew Vanderburg's keplersplinev2 to run: https://github.com/avanderburg/keplersplinev2.

## remake_pdf.ipynb

Remakes the posterior distribution functions (PDFs) from EXOFAST outputs for fits that are bimodal in stellar mass and age or for the radii of planets in grazing transit configurations.

## RV_mini_pipeline.py

A simple script that collects TRES or CHIRON radial velocity data from http://tess.exoplanets.dk/ and formats it for EXOFASTv2

## generate_priors_TOI.ipynb

A script that generates a spectroscopic metallicity prior, TESS dilution prior, and planetary ephemeris starting points for a given TESS Object of Interest.

## create_tess_table.ipynb

A short notebook to generate a TESS observations summary table for a given set of TIC and TOI ids. By default, it assumes that you are using SPOC data where available, and QLP when SPOC data are not available.

## create_followup_table.ipynb

A short notebook to generate a follow-up time-series observations summary table from information posted to ExoFOP.

## setup_fit.ipynb

A master notebook that can be used to set up an EXOFASTv2 fit from start to finish.
