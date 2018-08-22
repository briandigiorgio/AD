#!/usr/bin/env python

import time
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors, rc
from matplotlib.ticker import MultipleLocator, NullFormatter
import argparse
import scipy.stats as stats
from scipy.optimize import curve_fit
from adtools import *


desc = 'finds correlation between bulge size and deviation from Mi vs AD fit'

parser = argparse.ArgumentParser(description = desc)
parser.add_argument('-f', '--file', type = str, help = 'fits file to be read' +
                    'should be of the form SPX-GAU-MILESHC-composite_*.fits')
args = parser.parse_args()

if __name__ == '__main__':

    #assign variables
    print("Loading SPX file...")
    plate,ifu,Mi,Mie,Nmr,Nmre,gasrc,gasrce,ad,ade = read_ad_mag(args.file)

    ad_lim = [0,20000]
    bulge_lim = [0,1]
    Re = args.file[-11:-5]

    #load extra files
    print("Loading other files...")
    drpall = load_fits('drpall-v2_0_1.fits')
    bulgedata = load_fits('prob_table_nb4.fits')

    print("Matching numbers...")
    #creates array of prob_table galaxies' indices in drpall file
    drptoprob = match_cats(drpall['ifura'], drpall['ifudec'],
                    bulgedata['FITTED_RA'], bulgedata['FITTED_DEC'])

    #does the same thing for SPC-GAU... files to drpall file
    plateifuspx= np.asarray([plate[i] + ifu[i] for i in range(len(plate))])
    drpplate = drpall['plate'].astype(str)
    drpifu = drpall['ifudsgn'].astype(str)
    plateifudrp = np.asarray([drpplate[i]+drpifu[i] 
                  for i in range(len(drpplate))])

    spxtodrparray = [np.argmax(np.abs(plateifudrp == plateifuspx[i])) for i in 
               range(len(plateifuspx))]
    spxtodrp = np.asarray(spxtodrparray)

    [a,b], pcov = curve_fit(exponential, Mi, ad, sigma = ade)
    print([a,b])
    plt.semilogy(Mi, ad, 'b.')
    plt.errorbar(Mi, ad, yerr = ade, xerr = Mie, fmt = '.')
    t = np.linspace(-24,-17,100)
    plt.semilogy(t, exponential(t,a,b), 'r-')
    plt.show()

    dAD = ad - exponential(Mi,a,b)

    #sets limits
    axes = plt.gca()
    axes.set_xlim(bulge_lim)
    axes.set_ylim((-50000,50000))
    
    bulge = bulgedata['BULGE_TO_TOT_I'][drptoprob[spxtodrp]]
    bulgee = bulgedata['BULGE_TO_TOT_I_ERR'][drptoprob[spxtodrp]]
    t = np.linspace(0,1,100)
    [d,e], pcov2 = curve_fit(line, bulge, dAD, sigma = ade)
    plt.plot(t, e + d*t, 'r-')
    plt.errorbar(bulge, dAD,  yerr = ade, fmt = '.')
    plt.title("AD Deviation from Trend for Different Bulge Fraction")
    plt.xlabel("Bulge Fraction")
    plt.ylabel("AD Deviation from Fit (m/s)")
    #print([d,e])

    plt.savefig("dADbulge%s.png" % args.file[-11:-7])
    plt.show()   
