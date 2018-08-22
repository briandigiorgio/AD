#!/usr/bin/env python

import time
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors, rc
import argparse
import scipy.stats as stats
from scipy.optimize import curve_fit
import matplotlib.colors as colors
import matplotlib.cm as cm
from adtools import *

desc = 'finds correlation between Sersic index and deviation from Mi vs AD fit'

parser = argparse.ArgumentParser(description = desc)
parser.add_argument('--median', default = False, action = 'store_true')
parser.add_argument('--save', default = False, action = 'store_true')
args = parser.parse_args()

if __name__=='__main__':
    start = time.time()
    #loads all of the files from various radii into corresponding variables
    print("Loading files...")
    m25,ad25,ade25=read_mass('SPX-GAU-MILESHC-composite_0.25Re.fits')
    m50,ad50,ade50=read_mass('SPX-GAU-MILESHC-composite_0.50Re.fits')
    m75,ad75,ade75=read_mass('SPX-GAU-MILESHC-composite_0.75Re.fits')
    m100,ad100,ade100=read_mass('SPX-GAU-MILESHC-composite_1.00Re.fits')
    m125,ad125,ade125=read_mass('SPX-GAU-MILESHC-composite_1.25Re.fits')

    #combines variables into easier to use arrays
    marray = np.asarray([m25,m50,m75,m100,m125])
    adarray = np.asarray([ad25,ad50,ad75,ad100,ad125])
    adearray = np.asarray([ade25,ade50,ade75,ade100,ade125])
    xlabels = [0.25,0.5,0.75,1.00,1.25]

    print("Finding quintiles...")
    #finds quintiles for each of the data sets
    quintiles = np.asarray(list(map(lambda x: np.percentile(x, 
        [20,40,60,80,100]), [m25,m50,m75,m100,m125])))
    print(quintiles)

    #empty arrays to fill with data
    qarray = [[[] for i in range(5)] for i in range(5)]
    means = [[[] for i in range(5)] for i in range(5)]
    weights = [[[] for i in range(5)] for i in range(5)]
    stdevs = [[[] for i in range(5)] for i in range(5)]
    cmap = make_cmap(5, 'inferno')

    print("Making cuts, calculating, and plotting...")
    #loops through all of the quintiles of all of the radii
    #cuts the original data so the correct data is in the correct quintile
    for i in range(5):
        
        #making the cuts according t the quintiles
        q1cut = (marray[i]<quintiles[i][0])
        q2cut = (marray[i]>quintiles[i][0])*(marray[i]<quintiles[i][1])
        q3cut = (marray[i]>quintiles[i][1])*(marray[i]<quintiles[i][2])
        q4cut = (marray[i]>quintiles[i][2])*(marray[i]<quintiles[i][3])
        q5cut = (marray[i]>=quintiles[i][3])
        cuts = [q1cut,q2cut,q3cut,q4cut,q5cut]

        #takes mean of each quintile of each radius
        for j in range(5):

            #applying the appropriate cut to each quintile, calculate error
            qarray[i][j]=adarray[i][cuts[j]]
            weights[i][j] = np.power(adearray[i][cuts[j]], -2)
            means[i][j], stdevs[i][j] = weighted_avg_std(qarray[i][j],
                weights[i][j])
            means[i][j] = np.average(qarray[i][j])

            #finds user's desired quantity for each quintile
            if args.median:
                means[i][j] = np.median(qarray[i][j])
                #errorarray = 1.253*errorarray #coefficient from google search
                #stdevs[i] = 1.253*np.asarray(stdevs[i]) 
                plottype = 'Median'
            else:
                #meanarray[i][j] = np.mean(qarray[i][j])
                plottype = 'Mean'

        plt.errorbar(xlabels,means[i],yerr=stdevs[i], c = cmap[i],
                label="Mass Quintile %d"%(i+1))
        #plt.plot(xlabels,means[i], c = cmap[i], label="Mass Quintile %d"%(i+1))

    '''
    cmap = make_cmap(100, 'inferno')
    for j in range(100):
        plot = [adarray[i][j] for i in range(5)]
        plt.plot(xlabels, plot, c = cmap[j])
    '''

    plt.legend()
    axes = plt.gca()
    axes.set_xlim((.2,1.3))
    axes.set_ylim((0,120))
    plt.xticks(np.linspace(.25,1.25,5))
    plt.grid(True)
    plt.title(plottype + ' AD Profiles for Increasing Mass Quintiles')
    plt.xlabel('Re')
    plt.ylabel(plottype + ' Asymmetric Drift (m/s)')

    if args.save:
        plt.savefig("adquints.png")
    print("Time taken: %f seconds" % (time.time() - start))
    plt.show()
