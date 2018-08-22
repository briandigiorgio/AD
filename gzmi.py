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
    lier = loadfits('gal_list_v2_0_1_bpt_classify3.fits')
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
    bins = np.linspace(0,1,6)
    x = np.linspace(-23,-17,100)

    #more useful arrays
    c = make_cmap(len(bins), 'gnuplot')
    popts = np.zeros((len(Re), len(bins)))
    pcovs = np.zeros((len(Re), len(bins)))

    plt.figure(figsize=(8,12))
    for j in range(len(Re)):
        #get plate/ifu data for matching
        plate = spx[j]['plate']
        ifu = spx[j]['ifudesign']
        plateifuspx = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])

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
        ad = np.log10(spx[j]['ad2_em'][gzi])
        ade =np.log10(spx[j]['ad2_se'][gzi])
        Mi = spx[j]['elpetro_absmag'][:,5][gzi]
        Mie= spx[j]['elpetro_abmerr'][:,5][gzi]
        pcs = gz['P_CS_DEBIASED'][match['GZ1_INDX'][i[gzi]]]

        plt.subplot(321+j)
        for k in range(len(bins)-1):
            cut = (pcs > bins[k])*(pcs < bins[k+1])

            #plt.errorbar(Mi[cut], ad[cut], xerr=Mie[cut], yerr=ade[cut], 
            #        fmt='.', c=c[k], label = bins[k+1])#, alpha = .2)
            plt.plot(Mi[cut], ad[cut], '.', c=c[k], label = bins[k+1])
            popt,pcov = curve_fit(line, Mi[cut], ad[cut])#,
            #        sigma=ade[cut], maxfev= 10000)
            #print("%s: %d" % (k, np.sum(cut)))
            plt.plot(x, line(x, popt[0], popt[1]), c=c[k])
            popts[j,k] = popt[1]
            pcovs[j,k] = pcov[1,1]
            #plt.legend()
            plt.xlabel('Mi')
            plt.ylabel('AD')
            plt.title('%s Re' % Re[j])
            plt.grid(True)
            ax = plt.gca()
            #ax.set_yscale("log")
            ax.set_xlim((-23,-17))
            ax.set_ylim((2,5))

    for q in range(len(bins)-1):
        plt.subplot(326)
        r = np.asarray(Re).astype(float)
        plt.plot(r, popts[:,q], c=c[q], label = bins[q])
        #plt.fill_between(r, popts[:,q]-pcovs[:,q], popts[:,q]+pcovs[:,q],
        #        color=c[q], alpha = .2)
        plt.grid(True)
        plt.legend(loc = 0, fontsize='small')
        plt.xlabel('Re')
        plt.ylabel('Slope')
        plt.title('Slope for Different Radii')
    plt.tight_layout()
    plt.show()


    #not that eventful: groupings are relatively identically distributed
    #not that much variation in slope between them
