import astropy.units as u
import numpy as np
from astroquery.vizier import Vizier
import pandas as pd

import numpy as np 
import pandas as pd
import astropy.units as u
import requests

def grab_dilution(ticid):
    '''
    ticid: string containing the target's TIC ID. E.g. 'TIC 129004637' or '129004637'
    '''
    if ticid[:3] != 'TIC':
        ticid = 'TIC ' + ticid # adding 'TIC ' to the front of the ticid for Vizier to understand

    columns = ['TIC', 'Rcont']
    Vizier_search = Vizier(columns=columns)
    result = Vizier_search.query_object(ticid, catalog='IV/39/tic82', radius=2 * u.arcsec)[0]
    target_df = result.to_pandas()
    contam_ratio = target_df.Rcont[0]
    dilute_sigma = 0.1 * (contam_ratio / (1 + contam_ratio))
    print(f'Dilution prior: 0.0 +/- {np.round(dilute_sigma, 7)}')

    return dilute_sigma

def grab_metallicity(ticid, sigma=1, weighted=True):
    '''
    Generates a metallicity prior from the mean and standard deviation of spectroscopic metallicities available on ExoFOP.

    TICID: The object's TESS Input Catalog identifier. Ex: for TIC 81247740, enter '81247740'
    sigma: The number of standard deviations that you would like as your prior width.
    weighted: boolean to choose whether or not the average metallicity is weighted by the SNR of the spectrum
    '''

    url = 'https://exofop.ipac.caltech.edu/tess/target.php?id=' + ticid + '&json'
    page = requests.get(url)

    n_feh = page.json()['stellar_parameters'][4]['prov_num']

    feh = []
    snr = []
    for i in range(int(n_feh)):
        feh.append(float(page.json()['stellar_parameters'][5+i]['met']))
        snr.append(float(page.json()['stellar_parameters'][5+i]['snr']))
    feh = np.array(feh)
    snr = np.array(snr)

    if weighted==False:
        mean_feh = feh.mean()
    elif weighted==True:
        mean_feh = np.average(feh, weights=snr)
    width = sigma * np.std(feh)
    print(f'Spectroscopic [Fe/H] prior: {np.round(mean_feh, 4)} +/- {np.round(width,4)}')

    return mean_feh, width

def grab_all_priors(TOI, feh_sigma=1, feh_weighted=True):
    '''
    Grabs metallicity and dilution priors and ephemeris starting points for a given TOI and outputs it into a text file in EXOFASTv2 format.

    TOI: the TOI identifier for the target. E.g. 'TOI-3919'
    feh_sigma: the number of standard deviations to use as your metallicity prior width
    feh_weighted: whether or not to weight your average metallicity by the spectral SNR
    '''

    # Collect TICID and starting points
    url="https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"
    TOI_df=pd.read_csv(url, delimiter='|', index_col=1)
    TOI_id = TOI.replace('TOI-','') + '.01'
    ticid = 'TIC ' + str(TOI_df.loc[float(TOI_id)]['TIC ID'])

    depth = TOI_df.loc[float(TOI_id)]['Depth (ppm)']

    period = TOI_df.loc[float(TOI_id)]['Period (days)']
    tc = TOI_df.loc[float(TOI_id)]['Epoch (BJD)']
    rp_rstar = np.sqrt(depth * 10**(-6))

    print(f'Starting points:\nTc = {tc}\nPeriod = {period}\nRp/Rs = {rp_rstar}')
    feh, feh_width = grab_metallicity(ticid, sigma=feh_sigma, weighted=feh_weighted)
    dilute_sigma = grab_dilution(ticid)

    priorstring = f'# spectroscopic metallicity\nfeh {feh} {feh_width}\n# b\ntc_0 {tc}\nperiod_0 {period}\np_0 {rp_rstar}\ncosi_0 0.001\n# dilution\ndilute_0 0.0 {dilute_sigma}'

    with open(f'priors_{TOI}.txt', "w") as text_file:
        text_file.write(priorstring)