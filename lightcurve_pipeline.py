import numpy as np
import pandas as pd
import sys
import astropy
import astropy.units as u
#from keplersplinev2 import keplersplinev2
#from keplersplinev2 import *
from keplersplinev2 import choosekeplersplinev2
import lightkurve as lk
from lightkurve import search_targetpixelfile
import scipy as scipy
from scipy import stats
import matplotlib.pyplot as plt
import os
from astropy.time import Time

# only need these if wanting to look up values through code rather than entering them manually
import requests
from bs4 import BeautifulSoup as bsoup
import re # for parsing through ExoFOP (unconfirmed planets)

# general functions that can be used outside of making lcs

def bindata(x,y,binwidth):  
    
    '''
    This is a custom binning function. 
    x: time or phase array
    y: flux array
    binwidth: width of bins - I have NOT included astropy units here so make sure it's in the correct units!!!
    '''

    xbin = np.arange(np.nanmin(x),np.nanmax(x),binwidth)
    ybin = np.ones(len(xbin))
    #errbin = np.ones(len(xbin))
    
    for i in range(len(xbin)):
        mask = np.abs(x-xbin[i]) < binwidth/2
        #print(mask)
        ybin[i] = np.nanmedian(y[mask])

        #errbin[i] = np.sqrt(np.sum(err[mask]**2))
        
    # in case there are empty bins, exlcude these     
    nanmask = ~np.isnan(ybin)
    ybin = ybin[nanmask]
    xbin = xbin[nanmask]
    
    #print('nanmask: ',nanmask)
    return xbin,ybin#,errbin 

def calc_t14(Rp,Rstar,b,period,sma,inc):
    
    '''
    This returns the transit duration in hours, for use if it's not listed on NEA/EXOFOP etc.
    Astropy units ARE used here so make sure to include units on all inputs!

    
    Rp: planet radius
    Rstar: star radius
    b: impact parameter   --- if this isn't listed, set b = 0 to calculate longest possible transit duration
    Period should be self explanatory
    sma: semi-major axis
    inc: inclination (make sure to include deg or rad astropy units for input)
    '''
    
    C = (1.+(Rp/Rstar))**2 - b**2
    t14 = period/np.pi * np.arcsin(Rstar*C/sma/(np.sin(inc.to('rad'))).value)/u.rad
    
    return t14.to('hour')

def dilute(contam):
    
    '''
    Returns dilution term for a given TESS contamination ratio (listed at top left on EXOFOP).
    Include this in the prior file as: dilute_1 0.0 dilute_term
    '''
    
    return (contam/(1+contam))*0.1

# these functions are to make it easier to grab parameter values from NEA 

def find_num(number):
    
    '''
    This takes in a number and returns the main decimal number without the errors
     e.g. enter 3.23 +1.2-0.3, will return 3.23
    
    For use when grabbing values from NEA as they are in this format.
    '''
    
    ii=0    
    while ii < len(number):
        try:
            float(number[ii])
            
        except:
            if number[ii] != '.':
                break
        ii += 1                
    return float(number[:ii])


def get_array(arr):
    
    '''
    Appends the parameter entries that have values in the NEA table.
    
    '''
    
    array = []
    for ii in arr:
        thisnum = ii.text.strip()
        if thisnum != '---':
            array.append(thisnum)
    return array

def grabParam(paramName,result_list):
    
    '''
    Get parameter value from results list returned by scraping NEA page. 
    Returns the parameter with astropy units applied.
    '''

    results = result_list
    
    if paramName == 'Tc':
        title = "Time of Conjunction (Transit Midpoint)"
    elif paramName == 'Depth':
        title = 'Ratio of Planet Radius to Stellar Radius'
    else:
        title = paramName
    
    # this is grabbing the value for the specified parameter as long as it exists
    # pretty sure there's a much more efficient way of doing this
    
    if results.find('div',title=title) != None:
        findparam = results.find('div',title=title).parent.parent
        rows = findparam.find_all('td')
        param_list = get_array(rows)
        param = (find_num(param_list[0]))
    else:
        print('Parameter not found: ',title,'. Please specify parameter value and run manually.')
           
    if title == 'Orbital Period':
        param = param * u.day
    elif title == 'Semi-Major Axis':
        param = param * u.au
    elif title == 'Inclination':
        param = (param * u.deg).to('rad')
    elif title == 'Transit Duration':
        param = param * u.hour 
    elif paramName == 'Depth':
        param = param **2
    
    return param

def get_params(system,planet,published=True,target=None):
    
    '''
    Finds and returns all entries for the period, tc and transit duration on NEA for a given target.
    system: system name in K2-## format
    planet: letter of planet in system
    published: Boolean for whether or not the planet has been published. Reverts to ExoFOP if false.
    target: name of target system, e.g. 'TIC ###'
    '''
    if published == True:
        url = 'https://exoplanetarchive.ipac.caltech.edu/overview/' + system 
        page = requests.get(url)
        soup = bsoup(page.content,"html.parser")
    
        # add statement to get list of all planets and reformat as K2-##-#
        # maybe make another function to grab each column? unless there's an easier way
    
        # get planet data table and extract period, tc and duration from this
    
        results = soup.find(id = "planet_data_" + system + '-' + planet)
    
        # currently this makes an array of all entries for a given parameter, then grabs the first one which
        # is the most recent. might want to change this in future to grab param with smallest errors or something.
        # or maybe have a flag to print them out and have the user select?
    
        period = grabParam('Orbital Period',results)
        tc = grabParam('Tc',results)
        t14 = grabParam('Transit Duration',results)
        depth = grabParam('Depth',results)
        
    else: # added section for unconfirmed planets
        TIC = target.replace('TIC ', '')
        url = 'https://exofop.ipac.caltech.edu/tess/download_planet.php?id=' + TIC
        
        page = requests.get(url) # downloads raw text of the parameter table

        paramtable = str(page.content)
        
        # the following line uses a regular expression to grab the params. If ExoFOP changes how they format the table, this will
        # need to be changed. For a visual representation of regular expressions, see https://regex101.com/.
        parameters = re.findall(\
        "TOI \d+\.\d+\|+(\d*\.*\d*)\|\d*\.*\d+\|(\d*\.*\d+)\|\d*\.*\d*\|\d*\.*\d*\|\d*\.*\d*\|(\d*\.*\d*)\|\d*\.*\d*\|(\d*\.*\d*)", paramtable)
        if parameters == []:
            parameters = re.findall(\
            "\|+(\d*\.*\d*)\|\d*\.*\d+\|(\d*\.*\d+)\|\d*\.*\d*\|\d*\.*\d*\|\d*\.*\d*\|(\d*\.*\d*)\|\d*\.*\d*\|(\d*\.*\d*)\|\d*\.*\d+\|\d*\.*\d+\|\d*\.*\d+\|\d*\.*\d+\|\d*\.*\d+\|(\d*\.*\d+)\|", paramtable)
        paramlist = list(parameters[0]) # converts a list containing one tuple to a list of all 4 parameters
        
        period = float(paramlist[1]) * u.day # in JD
        tc = float(paramlist[0]) # in JD
        t14 = float(paramlist[3]) * u.hr # in hrs
        if paramlist[2] == '': # accounting for cases where only Rp/R* is reported
            Rp_Rstar = float(paramlist[4]) # Rp/R*
            depth = Rp_Rstar**2
        else:
            depth = float(paramlist[2]) / 1e6 # unitless (division by 1e6 is removing ppm)

        
    return period, tc, t14, depth

# main function for getting lightcurve

def create_light_curve(target, author, sector, period=None, duration=None, tc=None, targetname=None, 
    exposure = None,multisector = False, save = False, plot=True, binlc=False, binwidth = None, auto=False,
    system = None, planet = None, qualityFlag = False, depth = None, published = True, omit_transit_index=None,
    outputdir = 'None', outputfiles = 'all'):
    
    '''
    Create light curve files - currently for TESS only. If only a single sector is needed, use 
    create_light_curve directly. If lightcurves for all sectors/cadences are needed, use multi_sector.
    
    target: name of target system, e.g. 'TIC ###'
    author: pipeline name, e.g. 'SPOC'  
    sector: TESS sector
    period: planet period in days
    duration: transit duration in hours
    tc: time of transit midpoint as BJD
    targetname: optional name used for naming files. if not specified default to target
    exposure time: observing cadence. do not need to specify, used for multisector or if you want specific LC
    multisector: flag for creating lcs for multiple sectors. 
    save: flag to save lightcurves to text files    
    plot: flag to show plots in notebook
    binlc: flag to bin the lightcurve(s)
    binwidth: width of bins used if binning
    auto: automatically search the NEA for values for period, tc and transit duration
    system: system name in K2-## format. only needed if using automatic search for values, has to match NEA.
    planet: planet letter. only needed if using automatic NEA search for values.
    quality flag: if true only uses Lightkurve data that has zero bad quality flags. good for removing major spikes
    depth:
    published: Boolean for whether or not the planet has been published. Reverts to ExoFOP if false.
    omit_transit_index: the index of the transit(s) you wish to remove from the full lightcurve files
    outputdir: A string of the path to the directory where the files will be saved
    outputfiles: An array of the files that you wish to be generated. Options include: plot, fullnotflat, fullflat, fullnotrans, slimflat, and individual. Defaults to 'all'
    '''
    
    ############## 
    # added - erica
    
    # If running for multiple sectors, specify exposure time in case of multiple cadences in one sector
    # Can also specify exposure time for a single lc

    #Finding and downloading target information and light curve
    targetinfo = lk.search_lightcurve(target, author = author, sector = sector)    
       
    if multisector == True or (multisector == False and exposure != None):
        targetinfo = targetinfo[np.where(targetinfo.table['t_exptime']==exposure)].download()
    
    # if no exposure time is specified download first lk entry by default
    elif multisector == False and exposure == None:
        exposure =targetinfo.table['t_exptime'][0]
        if len(targetinfo) > 1:
            print('Warning! There are multiple entries for this sector. Only using first entry.')                
        targetinfo = targetinfo.download()
        
    targetinfo = targetinfo.remove_nans()
    
    if qualityFlag == True:
        # flag to only use data that has no bad data quality flags
        targetinfo = targetinfo[np.where(targetinfo.quality==0)]  
    
    # just checking if it's doing what i want
    print('Sector: %s. Exposure time: %s seconds.' %(sector, int(exposure)))
    
    if author == 'QLP'or 'TESS-SPOC':
       
        bjd = targetinfo.meta['BJDREFI']+targetinfo.meta['TSTART']
        date = Time(bjd,format='jd').fits[:10].replace('-','')
    
    else:
        date = targetinfo.meta['DATE'].replace('-','')
    
    # if target name isn't specified use input target by default
    if targetname == None:
        targetname = target
    
    if auto == True:
        period,tc, duration,depth = get_params(targetname,planet,published=published,target=target)
        
    
    #################
    
    #Determining Period, Transit Duration, and Tc
    period = period.to('day').value
    transit_duration = duration.to('hour').value
    tc = tc
    print('period: ',period,'duration: ',duration,'tc: ',tc) 
    #Creating a buffer for transit plotting
    buffer = 0.010 * transit_duration
    #print('buffer: ',buffer)
    
    #Creating target time and flux information
    time = np.array(targetinfo.time.value + 2457000)
    flux = np.array(targetinfo.flux.value)/np.median(np.array(targetinfo.flux.value))
    
    #Phase folding the light curve to identify transit properties
    phase = (time-tc)/period-np.floor((time-tc)/period)
    gt5 = phase > 0.5
    phase[gt5] = phase[gt5]-1.0
    
    #Determining the size of the phase transit
    transit_size_phase = ((transit_duration/24)/period+(buffer/24)/period)*3

    #Indicating where transits are located within the light curve
    in_transit = np.where((phase>-transit_size_phase/2) & (phase<transit_size_phase/2))
    out_of_transit = np.where(~((phase>-transit_size_phase/2) & (phase<transit_size_phase/2)))
    
    #Creating masks to exclude transits from the flattening process
    input_mask = np.ones_like(phase, dtype=bool)
    input_mask[in_transit] = False 
    
    #Flattening the light curve for any stellar variability
    flatten_masked, metadata = choosekeplersplinev2(time, flux, input_mask = input_mask, return_metadata = True)
    flat_flux = flux/flatten_masked

    
    binflag = 'unbinned'
    
    # If binning data, overwrite time, phase, flux and flat flux with binned versions
    if binwidth != None:
        print('Binning to :', binwidth)
        binflag = 'binned-%s' %str(binwidth).replace(' ','')
        binwidth = binwidth.to('day').value
        bintime, flux = bindata(time,flux,binwidth)
        bintime, flat_flux = bindata(time,flat_flux,binwidth)
        time = bintime
        phase = ((time)-tc)/period-np.floor(((time)-tc)/period)
        gt5 = phase > 0.5
        phase[gt5] = phase[gt5]-1.0
        
        in_transit = np.where((phase>-transit_size_phase/2) & (phase<transit_size_phase/2))
        out_of_transit = np.where(~((phase>-transit_size_phase/2) & (phase<transit_size_phase/2)))   
    
    in_transit_array = np.array(in_transit) # turning in_transit into an array from a one-element tuple

    #Separating the individual transits
    transit_times = []
    lc_transits = []
    array_names = []
    in_transit_n = []
    omit_transit_mintimeindex = []
    omit_transit_maxtimeindex = []
    floored_time = np.floor(((time)-tc)/period)
    unique_epochs = np.unique(floored_time)
    
    #print('unique_epochs: ',unique_epochs)
    
    for val in unique_epochs:                #Attaching the actual transit times to their location indicators
        t_time = tc + (val*period)
        transit_times.append(t_time)
    #print('transit_times: ', transit_times)
    for j in transit_times:                  #Only including transits within the observation timeframe
        if j > time[0] and j < time[-1]:
            lc_transits.append(j)
    #print('lc_transits: ', lc_transits)
    for current_time in lc_transits:                        
        index = lc_transits.index(current_time)         #Creating an array for naming each available transit
        if index != omit_transit_index: # only append new names if it isn't an omitted transit
            transit_name = 'in_transit_%s' % index
            array_names.append(transit_name)
        
            #Creating an array for where each transit can be found
            #in_transit_n.append(np.where((time>i-3/2*transit_duration-buffer) & (time<i+3/2*transit_duration+buffer)))  
            low_cut = (1.5*transit_duration-buffer)/24.
            hi_cut = (1.5*transit_duration+buffer)/24.
        
            in_transit_n.append(np.where((time>current_time-low_cut) & (time<current_time+hi_cut)))
        elif index == omit_transit_index:
            # saving the min and max indices of the transits to be omitted
            low_cut = (1.5*transit_duration-buffer)/24.
            hi_cut = (1.5*transit_duration+buffer)/24.

            omit_transit_mintimeindex.append(np.min(np.where((time>current_time-low_cut) & (time<current_time+hi_cut))))
            omit_transit_maxtimeindex.append(np.max(np.where((time>current_time-low_cut) & (time<current_time+hi_cut))))
    if omit_transit_index != None:
        for i in range(len(omit_transit_mintimeindex)):
            # cutting individual bad transits from the array of "in transit" indices
            in_transit_array = np.setdiff1d(in_transit_array, range(omit_transit_mintimeindex[i], omit_transit_maxtimeindex[i]))
    
    in_transit = in_transit_array.tolist()
    # the following if statement is to mute a warning. This may(?) have to be changed in the future and there might be a better way to do this
    if omit_transit_index == None:
        in_transit = tuple(in_transit)
    
    #print('array_names: ',array_names)
    #print('in_transit_n: ', in_transit_n)
    
    #Combining both the name array and location array into a comprehensive dictionary 
    # (a dataset for each transit)
    final_transits = dict(zip(array_names, in_transit_n))  
        
    # per-point error calculation - specific to TESS
    error_value = scipy.stats.median_abs_deviation(flat_flux[out_of_transit],nan_policy='omit')/0.68   #error per point
    print('per-point error value', error_value)
    errors = np.full((len(time), 1), error_value, dtype=float)
    
    # comparison of rms to transit depth
    if auto == True:
        rms = np.sqrt(np.sum((1.-flat_flux[out_of_transit])**2.)/len(flat_flux[out_of_transit]))
        print('Transith depth: %f, RMS: %.4f, Depth/RMS: %.3f' %(depth, rms,depth/rms))

    if plot == True:

        #Plotting the flattened, masked flux over the original flux to see the effects of cleaning the data
        fig1, ax1 = plt.subplots(figsize=(16, 8))
        plt.scatter(time, flux,label='Unflat Flux')
        plt.scatter(time, flat_flux, facecolors = 'tan', s = 100, alpha = 0.5, edgecolors='#000000',label='Flat Flux')
        plt.title("%s: Sector %s, %ss %s - Flattened, Masked Flux Vs. Original Flux" % \
                  (targetname,sector,int(exposure),binflag)) 
        plt.xlabel("Time - 2457000 [BTJD days]")
        plt.ylabel("Normalized Flux")
        plt.legend()

        #Plotting the light curve of only the transits
        fig2, ax2 = plt.subplots(figsize=(16, 8))
        plt.scatter(time[in_transit], flat_flux[in_transit])
        plt.title("%s: Sector %s, %ss %s - Individual Transit Flux" % (targetname,sector,int(exposure),binflag))
        plt.xlabel("Time - 2457000 [BTJD days]")
        plt.ylabel("Normalized Flux")
        #plt.ylim(0.992,1.01)
        #plt.xlim(2.4584e6+40.25,2.4584e6+41)

        #Plotting Individual Transits
        for key, value in final_transits.items():   
            newvalue = np.array(value)
            fig3, ax3 = plt.subplots(figsize=(16, 8))
            ax3.ticklabel_format(useOffset=False)
            plt.scatter((phase[newvalue[0,:]]*period), flat_flux[newvalue[0,:]])
            plt.title("%s: Sector %s, %ss %s - %s" % (targetname,sector,int(exposure),binflag, key))
            plt.xlabel("Time Since Transit Center (Days)")
            plt.ylabel("Normalized Flux")
            #print('transit times:',time[newvalue[0]])
            plt.show()
            plt.close()

        #Plotting the phase folded light curve to easily observe transit
        fig4, ax4 = plt.subplots(figsize=(16, 8))
        plt.scatter(phase[in_transit], flat_flux[in_transit])
        plt.title("%s: Sector %s, %ss %s - Phase Folded Light Curve" % (targetname,sector,int(exposure),binflag)) 
        plt.xlabel("Phase")
        plt.ylabel("Normalized Fluxf")
        #plt.xlim(min(phase[in_transit]),max(phase[in_transit]))
        #plt.ylim(0.99,1.01)
        
        

    # Saving data files of masked light curve, phase folded light curve, and light curve 
    # of the transits by themselves
    
    if save == True:

        if outputdir=='None':

            if os.path.exists('./lc_output') == False:
                os.mkdir('./lc_output')
            
            if ('plot' in outputfiles) or (outputfiles == 'all'):
                fig4.savefig('lc_output/%s_S%s_phase_folded_%ss%s.png' % (targetname,sector,int(exposure),binflag),dpi=400, bbox_inches="tight",format='png',facecolor='white')
            
            if ('fullnotflat' in outputfiles) or (outputfiles == 'all'):
                np.savetxt('lc_output/n%s.TESS.TESS.%sFullNotFlat.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time, flux, errors], delimiter=' ')
            if ('fullflat' in outputfiles) or (outputfiles == 'all'):
                np.savetxt('lc_output/n%s.TESS.TESS.%sFullFlat.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time, flat_flux, errors], delimiter=' ')
            if ('slimflat' in outputfiles) or (outputfiles == 'all'):
                np.savetxt('lc_output/n%s.TESS.TESS.%sSlimFlat.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time[in_transit], flat_flux[in_transit], errors[in_transit]], delimiter=' ')
            if ('fullnotrans' in outputfiles) or (outputfiles=='all'):
                np.savetxt('lc_output/n%s.TESS.TESS.%sFullNoTrans.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time[out_of_transit], flux[out_of_transit], errors[out_of_transit]], delimiter=' ')
            
            # save individual transits to separate files
            if ('individual' in outputfiles) or (outputfiles == 'all'):
                for key, value in final_transits.items():
                    newvalue = np.array(value)
                    np.savetxt('lc_output/n%s.TESS.TESS.%sSlimFlat.S%s.%ss%s.%s.dat' % (date,targetname, sector,int(exposure),binflag, key), \
                            np.c_[(phase[newvalue[0,:]]*period), flat_flux[newvalue[0,:]]], delimiter=' ')

        else:
            if os.path.exists(outputdir+'/lc_output') == False:
                os.mkdir(outputdir+'/lc_output')

            if ('plot' in outputfiles) or (outputfiles == 'all'):
                fig4.savefig(outputdir+'lc_output/%s_S%s_phase_folded_%ss%s.png' % (targetname,sector,int(exposure),binflag),dpi=400, bbox_inches="tight",format='png',facecolor='white')
            
            if ('fullnotflat' in outputfiles) or (outputfiles == 'all'):
                np.savetxt(outputdir+'lc_output/n%s.TESS.TESS.%sFullNotFlat.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time, flux, errors], delimiter=' ')
            if ('fullflat' in outputfiles) or (outputfiles == 'all'):
                np.savetxt(outputdir+'lc_output/n%s.TESS.TESS.%sFullFlat.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time, flat_flux, errors], delimiter=' ')
            if ('slimflat' in outputfiles) or (outputfiles == 'all'):
                np.savetxt(outputdir+'lc_output/n%s.TESS.TESS.%sSlimFlat.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time[in_transit], flat_flux[in_transit], errors[in_transit]], delimiter=' ')
            if ('fullnotrans' in outputfiles) or (outputfiles=='all'):
                np.savetxt(outputdir+'lc_output/n%s.TESS.TESS.%sFullNoTrans.S%s.%ss%s.dat' % (date,targetname, sector,int(exposure),binflag), \
                        np.c_[time[out_of_transit], flux[out_of_transit], errors[out_of_transit]], delimiter=' ')
            
            # save individual transits to separate files
            if ('individual' in outputfiles) or (outputfiles == 'all'):
                for key, value in final_transits.items():
                    newvalue = np.array(value)
                    np.savetxt(outputdir+'lc_output/n%s.TESS.TESS.%sSlimFlat.S%s.%ss%s.%s.dat' % (date,targetname, sector,int(exposure),binflag, key), \
                            np.c_[(phase[newvalue[0,:]]*period), flat_flux[newvalue[0,:]]], delimiter=' ')