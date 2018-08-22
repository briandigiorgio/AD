#!/usr/bin/env python

import time
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors, rc
import argparse
import scipy.stats as stats
from scipy.optimize import curve_fit
from adtools import *

desc = 'finds correlation between Sersic index and deviation from Mi vs AD fit'

parser = argparse.ArgumentParser(description = desc)
parser.add_argument('-f', '--file', type = str, help = 'fits file to be read' +
                    'should be of the form SPX-GAU-MILESHC-composite_*.fits')
args = parser.parse_args()

if __name__ == '__main__':

    #assign variables
    print("Loading SPC file...")
    plate,ifu,Mi,Mie,Nmr,Nmre,gasrc,gasrce,ad,ade = read_ad_mag(args.file)

    ad_lim = [0,50000]
    ser_lim = [0,10]

    #load extra files
    print("Loading other files...")
    drpall = load_fits('drpall-v2_0_1.fits')
    sersicdata = load_fits('prob_table_ser.fits')

    print("Matching numbers...")
    #creates array of prob_table galaxies' indices in drpall file
    drptoser,d2d = match_cats(drpall['ifura'], drpall['ifudec'],
                        sersicdata['FITTED_RA'], sersicdata['FITTED_DEC'])
    print(len(drptoser))

    plt.hist(d2d.value)
    plt.show()
    bad = []
    for i in range(len(d2d)):
        if d2d[i].is_within_bounds('10s',None):
            bad += [i]
    drptoser = np.delete(drptoser, bad)
    print(len(drptoser), len(bad))
 
    print(drpall['ifura'][0:20])
    print(sersicdata['FITTED_RA'][drptoser[0:20]])

    #does the same thing for SPC-GAU... files to drpall file
    plateifuspx = np.asarray([plate[i] + ifu[i] for i in range(len(plate))])
    drpplate = drpall['plate'].astype(str)
    drpifu = drpall['ifudsgn'].astype(str)
    plateifudrp = np.asarray([drpplate[i]+drpifu[i] 
                  for i in range(len(drpplate))])

    spxtodrparray = [np.argmax(np.abs(plateifudrp == plateifuspx[i])) for i in 
          range(len(plateifuspx))]
    spxtodrp = np.asarray(spxtodrparray)
    #plt.hist(sersicdata['SERSIC_INDEX'][drptoser], label = 'ser')
    #plt.hist(np.delete(drpall['nsa_sersic_n'][drpall['nsa_sersic_n']>0],bad), label = 'drpall')
    #plt.legend()
    #plt.show()

    plt.plot(np.delete(drpall['nsa_elpetro_absmag'][:,4],bad),
        sersicdata['TOTAL_MAG_R'][drptoser], 'b,')
    axes = plt.gca()
    #axes.set_xlim((0,10))
    plt.show()

    axes = plt.gca()
    axes.set_xlim((.5,10))
    ser = sersicdata['SERSIC_INDEX'][drptoser[spxtodrp]]
    sere = sersicdata['SERSIC_INDEX_ERR'][drptoser[spxtodrp]]

    plt.show()
    plt.loglog(ser, ad, 'b.')
    plt.errorbar(ser, ad, xerr = sere,  yerr = ade, fmt = '.')

    Re = args.file[-11:-5]
    plt.title('Asymmetric Drift for Different Sersic Indices at ' + Re)
    plt.xlabel('Sersic Index')
    plt.ylabel('Asymmetric Drift (m/s)')
    plt.grid(True)
    axes.xaxis.grid(True, which = 'minor')

    #plt.savefig("adser%s.png" % args.file[-11:-7])
    plt.show()   
