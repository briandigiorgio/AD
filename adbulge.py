#!/usr/bin/env python

import time
import numpy as np
from astropy.io import fits
from astropy import coordinates
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import colors, rc
from matplotlib.ticker import MultipleLocator, NullFormatter
import argparse
import scipy.stats as stats
from scipy.optimize import curve_fit
from adtools import *

desc = 'finds correlation between bulge size and AD for fixed Mi'

parser = argparse.ArgumentParser(description = desc)
parser.add_argument('-f', '--file', type = str, help = 'fits file to be read' +
                    'should be of the form SPX-GAU-MILESHC-composite_*.fits')
parser.add_argument('-l', '--lower', type = float, help = 'lower Mi cut')
parser.add_argument('-u', '--upper', type = float, help = 'upper Mi cut')
args = parser.parse_args()

if __name__ == '__main__':

    #assign variables
    print("Loading SPC file...")
    plate,ifu,Mi,Mie,Nmr,Nmre,gasrc,gasrce,ad,ade = read_ad_mag(args.file)

    ad_lim = [0,20]
    bulge_lim = [0,1]
    Re = args.file[-11:-5]

    #cut the magnitudes to the specified range
    print("Making cut...")
    Micut = (Mi < args.upper) * (Mi > args.lower)

    #load extra files
    print("Loading other files...")
    drpall = load_fits('drpall-v2_0_1.fits')
    bulgedata = load_fits('prob_table_nb4.fits')

    print("Matching numbers...")
    #creates array of prob_table galaxies' indices in drpall file
    drptoprob = match_cats(drpall['ifura'], drpall['ifudec'], 
                bulgedata['FITTED_RA'], bulgedata['FITTED_DEC'])

    #does the same thing for SPX-GAU... files to drpall file
    plateifuspx = np.asarray([plate[i] + ifu[i] for i in range(len(plate))])
    drpplate = drpall['plate'].astype(str)
    drpifu = drpall['ifudsgn'].astype(str)
    plateifudrp = np.asarray([drpplate[i]+drpifu[i] 
                  for i in range(len(drpplate))])

    print(np.asarray([np.argwhere(plateifudrp == plateifuspx[i]) for i in
        range(len(plateifuspx))])[:,0,0])

    spxtodrparray = np.asarray([np.argwhere(plateifudrp == plateifuspx[i]) 
        for i in range(len(plateifuspx))])[:,0,0]
    spxtodrp = spxtodrparray[Micut]
    
    plt.plot(Mi[Micut], drpall['nsa_elpetro_absmag'][:,5][spxtodrp], 'b.')
    plt.show()

    #in case you are only working with a set Re, you can save/load this too
    #time saved is negligible
    ##np.save('spxtodrp', spxtodrparray)
    ##spxtodrp = np.load('spxtodrp.npy')[Micut]

    #cuts ad data and gets the correct bulge data through chain of indices
    #x and y are backwards
    y = ad[Micut]/1000
    x = bulgedata['BULGE_TO_TOT_I'][drptoprob[spxtodrp]]
    yerr = ade[Micut]/1000
    xerr = bulgedata['BULGE_TO_TOT_I_ERR'][drptoprob[spxtodrp]]

    print('Fitting...')
    for i in range(len(xerr)):
        if not yerr[i]:
            yerr[i] = .0001
    
    rho, p = stats.spearmanr(x,y)
    print(rho,p)
    #popt, pcov = curve_fit(line, x, y, sigma = yerr, absolute_sigma = True) 
    #print(popt, pcov)
    #points = np.linspace(0,20000,200)
    #plt.plot(points, points * popt[0] + popt[1], 'r-')

    #sets limits
    axes = plt.gca()
    axes.set_ylim(ad_lim)
    axes.set_xlim(bulge_lim)

    #plots
    print("Plotting...")
    plt.errorbar(x,y,xerr = xerr,yerr = yerr, fmt = '.')
    plt.grid(True)

    plt.title('Asymmetric Drift for Different Bulge Fraction at ' + Re
              +'\n for Galaxies with %d < Mi <%d' % (args.lower, args.upper))
    plt.xlabel('Bulge Fraction')
    plt.ylabel('Asymmetric Drift (km/s)')
    plt.savefig("adbulge%s%d%d.png"%(args.file[-11:-7], args.lower, args.upper))
    plt.show()
