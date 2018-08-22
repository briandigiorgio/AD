#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *
from astropy.coordinates import Angle
import time
from scipy.optimize import curve_fit

if __name__ == '__main__':
    #load all relevant files
    start = time.time()
    gz = loadfits('zoo2MainSpecz.fits.gz')
    match = loadfits('manga_catalog_match.fits.gz', i=2)
    spx25 = loadfits('SPX-GAU-MILESHC-composite_0.25Re.fits')
    spx50 = loadfits('SPX-GAU-MILESHC-composite_0.50Re.fits')
    spx75 = loadfits('SPX-GAU-MILESHC-composite_0.75Re.fits')
    spx10 = loadfits('SPX-GAU-MILESHC-composite_1.00Re.fits')
    spx12 = loadfits('SPX-GAU-MILESHC-composite_1.25Re.fits')
    col = 't03_bar_a06_bar_debiased'
    sp = 't04_spiral_a08_spiral_debiased'

    #make useful arrays
    spx = np.asarray([spx25, spx50, spx75, spx10, spx12])
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    popt = [[] for i in range(5)]
    pcov = [[] for i in range(5)]
    x = np.linspace(0,1,100)
    bins = np.linspace(0,1,6)
    c = make_cmap(len(bins-1), 'inferno')
    data = np.zeros((len(Re), len(bins)-1))
    errors = np.zeros_like(data)
    
    for j in range(len(Re)):
        #assign variables to relevant data
        plate = np.asarray(spx[j]['plate'])
        ifu = np.asarray(spx[j]['ifudesign'])

        #generates array matching galaxyzoo to drpall, adapted from kyle readme
        #does same function as python code below but many orders of mag faster
        i = np.asarray([np.arange(len(match['PLATE']))[(match['PLATE']==p) & 
                      (match['IFUDESIGN']==i)][0] for p,i in zip(plate,ifu)])
        '''
        index = []
        for i in range(len(ifu)):
            for j in range(len(match['PLATE'])):
                if match['PLATE'][j] == int(plate[i]): 
                    if match['IFUDESIGN'][j] == int(ifu[i]):
                        index += [j]
        i = np.asarray(index)
        np.save('spxtomatch', i)
        #i = np.load('spxtomatch.npy')
        '''

        #filter out bad values and pick out correct data
        gzi = (match['GZ2_INDX'][i] > -1)
        spirals = (gz[sp][match['GZ2_INDX'][i[gzi]]]>0.5)
        ad = spx[j]['ad2_em'][gzi][spirals]
        ade =spx[j]['ad2_se'][gzi][spirals]

        pcs = gz[col][match['GZ2_INDX'][i[gzi]]][spirals]
        #lt.errorbar(pcs,ad,yerr=ade,fmt='.',c=c[j],label='%s Re' % Re[j])
        popt[j], pcov[j] = curve_fit(line, pcs, ad) 
        #plt.semilogy(x, exponential(x, popt[j][0], popt[j][1]))

        for w in range(len(bins)-1):
            cut = (pcs > bins[w]) * (pcs < bins[w+1])
            data[j,w], errors[j,w] = weighted_avg_std(ad[cut], ade[cut])

        #plt.errorbar(bins[:len(bins)-1], means, yerr = stds, c=c[j], label = Re[j])
    for k in range(data.shape[1]):
        plt.errorbar(np.asarray(Re).astype(float), data[:,k], yerr=errors[:,k], 
                label = bins[k], c=c[k])
    plt.title('Galaxy Zoo Bar Likelihood vs. Mean AD for Different Re')
    plt.xlabel('Re')
    plt.ylabel('Asymmetric Drift')
    plt.legend(title = 'Bar Prob')
    plt.grid(True)
    print(popt)
    print("Time taken: ", time.time()-start)
    ax = plt.gca()
    ax.set_yscale("log")
    plt.show()

#baradsp: looking at mean ad magnitude for different bar probabilities for
#   galaxies with P(sp) > .5 yields an essentially flat curve so no variation
#   in ad due to bar
#
