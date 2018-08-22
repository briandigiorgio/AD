#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *
from astropy.coordinates import Angle
import time
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

if __name__ == '__main__':
    #load all relevant files
    start = time.time()
    lier = loadfits('gal_list_v2_0_1_bpt_classify3.fits')
    match = loadfits('manga_catalog_match.fits.gz', i=2)
    spx25 = loadfits('SPX-GAU-MILESHC-composite_0.25Re.fits')
    spx50 = loadfits('SPX-GAU-MILESHC-composite_0.50Re.fits')
    spx75 = loadfits('SPX-GAU-MILESHC-composite_0.75Re.fits')
    spx10 = loadfits('SPX-GAU-MILESHC-composite_1.00Re.fits')
    spx12 = loadfits('SPX-GAU-MILESHC-composite_1.25Re.fits')
    adr = False

    #make useful arrays
    spx = np.asarray([spx25, spx50, spx75, spx10, spx12])
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified')
    x = np.linspace(-23,-17,100)

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = np.nanmax(bpt) + 1
    bpt[np.isnan(bpt)] = len(types) - 1
    maxbpt = int(np.max(bpt))

    #more useful arrays
    c = make_cmap(int(np.max(bpt) + 1), 'gnuplot')
    popts = np.zeros((len(Re), maxbpt+1))
    pcovs = np.zeros((len(Re), maxbpt+1))

    plt.figure(figsize=(8,12))
    for j in range(len(Re)):
        #plt.figure(figsize=(8,12))
        #get plate/ifu data for matching
        plate = spx[j]['plate'].astype(str)
        ifu = spx[j]['ifudesign'].astype(str)
        plateifuspx = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])

        lierplate = lier['PLATE'].astype(str)
        lierifu = lier['IFUDESIGN'].astype(str)
        plateifulier = np.asarray([lierplate[i] + lierifu[i] 
            for i in range(len(lierplate))])

        #match the catalogs
        #for some reason the numpy version doesn't work but the python one does
        #spxtolier = np.asarray([np.argmax(plateifulier==plateifuspx[i])
        #                for i in range(len(plateifuspx))])

        spxtolier = np.zeros(len(plateifuspx))
        for l in range(len(plateifuspx)):
            for m in range(len(plateifulier)):
                if plateifuspx[l] == plateifulier[m]:
                    spxtolier[l] = m
        spxtolier = spxtolier.astype(int)

        #filter out bad values and pick out correct data
        ad = spx[j]['ad2_em']
        ade =spx[j]['ad2_se']
        harc = spx[j]['harc_em']
        harce = spx[j]['harc_se']
        Mi = spx[j]['elpetro_absmag'][:,5]
        Mie= spx[j]['elpetro_abmerr'][:,5]
        bad = np.where(np.isnan(np.log(harc*harce*ad*ade)))

        ad = np.delete(ad, bad)
        ade = np.delete(ade, bad)
        harc = np.delete(harc, bad)
        harce = np.delete(harce, bad)
        Mi = np.delete(Mi, bad)
        Mie = np.delete(Mie, bad)

        #harc[np.isnan(harc) or not harc] = 1
        if adr:
            ad = ad/harc
            ade = np.sqrt((ade/harc)**2 + ((ad*harce)/(harc**2))**2)
        fade = ade/ad
        ade = ade/ad
        ad = np.log10(ad)
        #ade = ad * fade


        plt.subplot(321+j)
        #for k in range(maxbpt+1):
        for k in [1,2]:
            cut = np.delete((bpt==k)[spxtolier], bad)
            #plt.subplot(321+k)
            plt.errorbar(Mi[cut], ad[cut], xerr=Mie[cut], yerr=ade[cut], 
                    fmt='.', c=c[k], label = types[k], alpha = .2)
            popt,pcov = curve_fit(line, Mi[cut], ad[cut],
                    sigma=ade[cut], maxfev= 10000)
            #print("%s: %d" % (k, np.sum(cut)))
            plt.plot(x, line(x, popt[0], popt[1]), c=c[k])
            popts[j,k] = popt[0]
            pcovs[j,k] = pcov[0,0]
            #plt.legend()
            plt.xlabel(r'$M_i$')
            plt.ylabel(r'$AD^2$')
            ax = plt.gca()
            ax.set_ylim((0,6))
            if adr:
                plt.ylabel(r'$AD^2/H_{rot}$')
                ax.set_ylim((-1,3))
            plt.title(r'%s $R_e$' % Re[j])
            plt.grid(True)
            #ax.set_yscale("log")
        plt.tight_layout()
        #plt.show()

    #for q in range(maxbpt+1):
    for q in [1,2]:
        plt.subplot(326)
        r = np.asarray(Re).astype(float)
        plt.plot(r, popts[:, q], c=c[q], label = types[q])
        plt.fill_between(r, popts[:,q]-pcovs[:,q], popts[:,q]+pcovs[:,q],
                color=c[q], alpha = .2)
        plt.grid(True)
        plt.legend(loc = 0, fontsize='small')
        plt.xlabel(r'$R_e$')
        plt.ylabel('Slope')
        plt.title('Slope for Different Radii')
        ax = plt.gca()
        ax.set_ylim((-.5,.1))
    plt.tight_layout()
    plt.show()
