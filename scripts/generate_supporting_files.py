import numpy as np
import pandas as pd

def generate_pro_file(toi, ticid, exptimes, tess_lcs, nplanets=1, nstars=1, outpath='./'):
    '''
    Generates the EXOFASTv2 procedures file to start a fit, assuming that transits and RVs are both being fit and that the RV fit is eccentric and has a linear slope.

    toi: TOI id (numbers only)
    ticid: TIC id (numbers only)
    exptimes: An array containing the exposure time in minutes of each lightcurve by index. E.g. [30, 10, 2, 2, 2]
    tess_lcs: An array containing booleans to determine which lightcurve is a tess lightcurve by index. E.g. [1, 1, 0, 0, 1]
    nplanets: The number of planets in the system.
    nstars: The number of stars in the system.
    outpath: Directory to send the generated procedures file.
    '''

    # Create ninterps array using Jason's recommended strategy
    ninterps = []
    for exptime in exptimes:
        if exptime == 30:
            ninterps.append(10)
        elif exptime == 10:
            ninterps.append(4)
        elif exptime < 10:
            ninterps.append(1)
        else:
            ninterps.append(1)
            raise Warning('Unexpected exposure time. Please manually edit the ninterp list and check the exptimes list.')
    
    exptimes_str = str(exptimes)
    ninterps_str = str(ninterps)
    tess_lcs_str = str(tess_lcs)


    boilerplate = f'''pro fittoi{toi}, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
nthread=nthread, outpath=outpath, maxtime=maxtime

exofastv2, nplanets={nplanets}, nstars={nstars}, tranpath='n2*.dat', $
        rvpath='*.rv',priorfile='{ticid}.priors', $
        prefix='./fitresults/{ticid}.', maxsteps=maxsteps, $
        nthin=nthin, circular=[0], fitrv=[1], fittran=[1], maxtime=maxtime, $
        debug=debug, verbose=verbose, nthread=nthread, mistsedfile='{ticid}.sed', /fitslope, $
    exptime={exptimes_str}, $
    ninterp={ninterps_str}, $
    fitdilute={tess_lcs_str}


end'''
    
    filename = f'fittoi{toi}.pro'
    with open(outpath + filename, 'w') as file:
        file.write(boilerplate)

def generate_SLURM_file(toi, ncpu = 20, runtime = 2, ttbuffer = 0.25, email = 'jschulte@msu.edu', outpath='./'):
    
    runtime_minutes = runtime * 24 * 60
    runtime_seconds = runtime_minutes * 60

    ttbuffer_seconds = ttbuffer * 24 * 60 * 60
    runtime_seconds_ttbuffer = int(runtime_seconds - ttbuffer_seconds)

    boilerplate = f'''#!/bin/bash --login

#SBATCH --account=exoplanet_lab # Set this to charge cpu-hours to buy-in account
# #SBATCH --nodelist="acm-007"    # Set this job to run on the buy-in node
#SBATCH --time={int(runtime_minutes)}              # Set the runtime limit for this (minutes)
#SBATCH --nodes=1               # Set to one node (for IDL license reasons)
#SBATCH --ntasks=1             # Number of tasks this job will run
#SBATCH --cpus-per-task={ncpu}       # Number of cores needed for each task
#SBATCH --mem-per-cpu=1G        # Memory required for each core
#SBATCH --mail-type=ALL         # SLURM automated messages to email
#SBATCH --mail-user={email}
#SBATCH --job-name=toi-{toi}     # Set a job name to easily track status
#SBATCH --open-mode=append.      # Set the output file to append mode
#SBATCH --output=FITLOG_%j.out # Standard output goes to this file

# Load the IDL module
module purge
module load IDL

# Setup the EXOFASTv2 environment
source /mnt/research/Exoplanet_Lab/source.sh

# run the fit

/opt/software-current/2023.06/x86_64/generic/software/IDL/idl88/bin/idl -e "fittoi{toi}, nthread={ncpu}, maxsteps=15000, nthin=700, maxtime={runtime_seconds_ttbuffer}"

# {runtime} days for SLURM, {runtime - ttbuffer} days for EXOFASTv2

wait

cd fitresults
makepdfs.sh
'''

    filename = f'SLURM_toi{toi}.sb'
    with open(outpath + filename, 'w') as file:
        file.write(boilerplate)