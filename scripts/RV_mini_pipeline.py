import pandas as pd
import numpy as np
import astropy.units as u
import requests


def download_TRES_data(ticid, username, password, output_filename='tres_rvs.txt', verbose=False):
    '''
    Download TRES or CHIRON data from Lars Buchhave's TRES website. It seems like this will download multi-orders only if they exist, but please check to be sure.

    ticid: string containing the TIC identifier of the target.
    username: your username to access the TRES website.
    password: your password to access the TRES website.
    output_filename: the filename that the output RV file should have.
    verbose: whether or not to print when login/download are successful.
    '''

    # URL of the login page
    login_url = 'http://tess.exoplanets.dk/Login.aspx'

    # convert ticid to T0 number
    T0_id = 'T0' + ticid.zfill(9)
    # URL of the text file you want to download
    file_url = f'http://exoplanets.dk/uploads/temp/{T0_id}.RV.txt'

    # Your login credentials
    payload = {
        'username': username,
        'password': password
    }

    # Create a session to persist the login state
    session = requests.Session()

    # Perform the login
    response = session.post(login_url, data=payload)

    # Check if the login was successful
    if response.status_code == 200:
        if verbose==True:
            print("Login successful!")
        
        # Download the file
        file_response = session.get(file_url)
        
        if file_response.status_code == 200:
            # Save the file locally
            with open(output_filename, 'wb') as file:
                file.write(file_response.content)
            if verbose==True:
                print("File downloaded successfully!")
        else:
            print("Failed to download the file.")
    else:
        print("Login failed.")

def remake_rvfile(path, units, output_filename='output.rv', verbose=False, download_rvs = False, ticid=None, username=None, password=None):
    '''
    Create a streamlined RV file from a TRES or CHIRON data file generated by Lars Buchhave's TESS website.

    path: string containing the path to the source data file.
    units: astropy units of the input RVs.
    output_filename: string containing the name of the processed RV file. This can contain a path to a different directory.
    verbose: a boolean for whether or not to print the processed RV file.
    download_rvs: a boolean for whether or not the RVs should be downloaded from Lars Buchhave's TRES/CHIRON website first.
    ticid: if download_rvs is True, the TIC identifier of the target must be provided.
    username: if download_rvs is True, your TRES/CHIRON site username must be provided.
    password: if download_rvs is True, your TRES/CHIRON site password must be provided.
    '''
    
    if download_rvs == True:
        download_TRES_data(ticid=ticid, username=username, password=password, output_filename=path, verbose=verbose)

    input_data = pd.read_csv(path, sep='\s+', header=0)

    vrad_units = np.array(input_data.vrad) * units
    svrad_units = np.array(input_data.svrad) * units
    
    data_new = pd.DataFrame({'Time (BJD)': input_data.BJD_UTC, 'RV (m/s)': vrad_units.to(u.m/u.s), 'RV error (m/s)': svrad_units.to(u.m/u.s)})

    with open(output_filename, 'w') as file:
        file.write('# Time (BJD), RV (m/s), RV error (m/s)\n')
        data_new.to_csv(file, sep=' ', index=False, header=False)

    if verbose==True:
        print('Processed file:\n', data_new)