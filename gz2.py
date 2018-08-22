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
    gz = loadfits('GalaxyZoo1_DR_table2.fits')
    match = loadfits('manga_catalog_match.fits.gz', i=2)
    spx25 = loadfits('SPX-GAU-MILESHC-composite_0.25Re.fits')
    spx50 = loadfits('SPX-GAU-MILESHC-composite_0.50Re.fits')
    spx75 = loadfits('SPX-GAU-MILESHC-composite_0.75Re.fits')
    spx10 = loadfits('SPX-GAU-MILESHC-composite_1.00Re.fits')
    spx12 = loadfits('SPX-GAU-MILESHC-composite_1.25Re.fits')

    #make useful arrays
    spx = np.asarray([spx25, spx50, spx75, spx10, spx12])
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    popt = [[] for i in range(5)]
    pcov = [[] for i in range(5)]
    x = np.linspace(0,1,100)
    c = make_cmap(len(Re), 'inferno')

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
        gzi = match['GZ1_INDX'][i] > -1
        ad = spx[j]['ad2_em'][gzi]
        ade =spx[j]['ad2_se'][gzi]

        pcs = gz['P_EL_DEBIASED'][match['GZ1_INDX'][i[gzi]]]
        #plt.errorbar(pcs,ad,yerr=ade,fmt='.',c=c[j],label='%s Re' % Re[j])
        popt[j], pcov[j] = curve_fit(line, pcs, ad) 
        #plt.semilogy(x, exponential(x, popt[j][0], popt[j][1]))

        bins = np.linspace(0,1,11)
        cuts = [[] for i in range(len(bins)-1)]
        means = np.zeros(len(cuts))
        stds = np.zeros(len(cuts))
        for i in range(len(cuts)):
            cuts[i] = (pcs > bins[i]) * (pcs < bins[i+1])
            means[i], stds[i] = weighted_avg_std(ad[cuts[i]], ade[cuts[i]])

        plt.errorbar(bins[:10], means, yerr = stds, c=c[j], label = Re[j])
    plt.title('Galaxy Zoo Elliptical Likelihood vs. Mean AD for Different Re')
    plt.xlabel('Debiased Spiral Probability')
    plt.ylabel('Asymmetric Drift')
    plt.legend()
    plt.grid(True)
    print(popt)
    print("Time taken: ", time.time()-start)
    ax = plt.gca()
    ax.set_yscale("log")
    plt.show()

#gzadcs.png: scatter plot shows slight downward trend
#gzadcsm.png: means confirm downward trend, show no meaningful diff between Re
#gzadelm.png: elliptical shows opposite
#shouldn't the more obviously spiral galaxies have higher AD?
#or maybe larger spirals look more elliptical
