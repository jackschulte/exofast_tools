pro fit[SYSTEM_NAME], debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
                nthread=nthread, outpath=outpath

  maxsteps=[MAXSTEPS]
  nthin=[NTHIN]

;; Fit using the Torres relation
exofastv2, nplanets=1, tranpath='n2*.dat', $
           rvpath='[RV_FILENAME]',priorfile='[PRIOR_FILENAME]', $
           prefix='./[FIT_DIR_NAME]/[SYSTEM_NAME].', maxsteps=maxsteps, $
           nthin=nthin, circular=[0], fitrv=[1],fittran=[1], $
           debug=debug, verbose=verbose, nthread=[NUM_THREADS], mistsedfile='[SED_FILENAME]', $
           /fitslope, exptime=[EXPTIME_LIST_MINUTES], ninterp=[NUM_INTERPS]


end
