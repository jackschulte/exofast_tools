pro fittoi4138, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
                nthread=nthread, outpath=outpath

  maxsteps=10000
  nthin=400

;; Fit using the Torres relation
exofastv2, nplanets=1, tranpath='n2*.dat', $
           rvpath='toi4138.TRES.rv',priorfile='toi4138.priors', $
           prefix='./fitresults_init/TOI-4138.MIST.SED.', maxsteps=maxsteps, $
           nthin=nthin, circular=[0], fitrv=[1],fittran=[1], $
           debug=debug, verbose=verbose, nthread=6, mistsedfile='toi4138.sed', $
           /fitslope, exptime=[30, 30, 30, 30, 30, 30, 30, 2, 2, 2, 2, 2, 2, 1, 1], ninterp=[10, 10, 10, 10, 10, 10, 10, 1, 1, 1, 1, 1, 1, 1, 1]


end
