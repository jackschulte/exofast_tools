import astropy.units as u
import numpy as np
from astroquery.vizier import Vizier
import pandas as pd

import numpy as np 
import pandas as pd
import astropy.units as u
import requests
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
import time
from io import StringIO
from selenium.webdriver.chrome.options import Options

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

def grab_metallicity(ticid, sigma=1, weighted=True, source='exofop', tres_username=None, tres_password=None, toi_id = None, verbose=False):
    '''
    Generates a metallicity prior from the mean and standard deviation of spectroscopic metallicities available on ExoFOP.

    TICID: The object's TESS Input Catalog identifier. Ex: for TIC 81247740, enter '81247740'
    sigma: The number of standard deviations that you would like as your prior width.
    weighted: Boolean to choose whether or not the average metallicity is weighted by the SNR of the spectrum.
    source: Either 'exofop' or 'tres'. 'exofop' pulls from the available metallicities on ExoFOP-TESS, while 'tres' pulls from Lars Buchhave's planet candidate site.
    tres_username: Your username for tess.exoplanets.dk
    tres_password: Your password for tess.exoplanets.dk
    toi_id: The target's TOI id. Ex: 'TOI-1855'
    verbose: A boolean to determine whether or not the program opens up the Chrome test browser to display what it's doing

    '''
    if source == 'exofop':
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
    elif source == 'tres':
        # Configure WebDriver (e.g., Chrome)
        chrome_options = Options()
        chrome_options.add_argument("--headless")  # Enable headless mode
        chrome_options.add_argument("--no-sandbox")  # Bypass OS security model (useful for Docker)
        chrome_options.add_argument("--disable-dev-shm-usage")  # Overcome limited resource problems

        if verbose==True:
            driver = webdriver.Chrome()
        else:
            driver = webdriver.Chrome(options=chrome_options)
        

        # Open the login page
        driver.get("http://tess.exoplanets.dk/Login.aspx")

        # Find and fill out the login fields
        username = driver.find_element(By.NAME, "ctl00$cpMainContent$tbUserName")
        password = driver.find_element(By.NAME, "ctl00$cpMainContent$tbPassword")

        username.send_keys(tres_username)
        password.send_keys(tres_password)
        password.send_keys(Keys.RETURN)

        # Wait for login to complete and navigate to the target page
        time.sleep(5)

        # Navigate to the page with the search field
        driver.get("http://tess.exoplanets.dk/Default.aspx")

        # Find the search field and perform the search
        search_field = driver.find_element(By.NAME, "ctl00$cpMainContent$tbSearch")
        search_field.send_keys(toi_id)
        search_field.send_keys(Keys.RETURN)

        # Wait for the search results to load
        time.sleep(5)

        # Find the link with href containing "Candidate_Edit" (the page for an individual target)
        links = driver.find_elements(By.TAG_NAME, "a")
        candidate_edit_link = None
        for link in links:
            href = link.get_attribute("href")
            if href and "Candidate_Edit" in href:
                candidate_edit_link = href
                break

        # Check if we found the correct link and navigate to it
        if candidate_edit_link:
            driver.get(candidate_edit_link)
        else:
            print("Candidate_Edit link not found")
            driver.quit()
            exit()

        # Wait for the page to load
        time.sleep(5)

        spc_table_button = driver.find_element(By.NAME, 'ctl00$cpMainContent$Button12')
        spc_table_button.click() # this opens a second tab

        # Wait for the page to load
        time.sleep(5)

        tabs = driver.window_handles 
        driver.switch_to.window(tabs[1])

        # You can extract text from a specific element or the entire page
        spc_table = driver.find_element(By.TAG_NAME, 'body').text

        # Create a DataFrame
        data = StringIO(spc_table)
        df = pd.read_csv(data, header=0, sep='\s+')
        df.dropna(inplace=True)
        df_spc = df[df.method=='SPC2.9']
        
        if weighted==False:
            mean_feh = df_spc.mh.mean()
        elif weighted==True:
            mean_feh = np.average(df_spc.mh, weights=df_spc.SNRe)
        width = sigma * np.std(df_spc.mh)

        # Quit the WebDriver
        driver.quit()
    print(f'Spectroscopic [Fe/H] prior: {np.round(mean_feh, 4)} +/- {np.round(width,4)}')

    return mean_feh, width

def grab_all_priors(TOI, feh_sigma=1, feh_weighted=True, outpath='.', source='exofop', tres_username=None, tres_password=None, verbose=False):
    '''
    Grabs metallicity and dilution priors and ephemeris starting points for a given TOI and outputs it into a text file in EXOFASTv2 format.

    TOI: the TOI identifier for the target. E.g. 'TOI-3919'
    feh_sigma: the number of standard deviations to use as your metallicity prior width
    feh_weighted: whether or not to weight your average metallicity by the spectral SNR
    outpath: the path to the generated prior text file. Defaults to the current working directory.
    source: Either 'exofop' or 'tres'. 'exofop' pulls from the available metallicities on ExoFOP-TESS, while 'tres' pulls from Lars Buchhave's planet candidate site.
    tres_username: Your username for tess.exoplanets.dk
    tres_password: Your password for tess.exoplanets.dk
    toi_id: The target's TOI id. Ex: 'TOI-1855'
    verbose: A boolean to determine whether or not the program opens up the Chrome test browser to display what it's doing
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
    feh, feh_width = grab_metallicity(ticid, sigma=feh_sigma, weighted=feh_weighted, source=source, tres_username=tres_username, tres_password=tres_password, toi_id = TOI, verbose=verbose)
    dilute_sigma = grab_dilution(ticid)

    priorstring = f'# spectroscopic metallicity\nfeh {feh} {feh_width}\n# b\ntc_0 {tc}\nperiod_0 {period}\np_0 {rp_rstar}\ncosi_0 0.001\n# dilution\ndilute_0 0.0 {dilute_sigma}'

    if outpath[-1] == '/':
        outpath = outpath[:-1]
    with open(f'{outpath}/priors_{TOI}.txt', "w") as text_file:
        text_file.write(priorstring)