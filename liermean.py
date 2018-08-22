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
    adr = True #set to normalize by gas rotation speed

    #make useful arrays
    spx = np.asarray([spx25, spx50, spx75, spx10, spx12])
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified','','','SF Control')
    x = np.linspace(-23,-17,100)

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    maxbpt = int(np.max(bpt))

    #parameters of interest
    interested = [8,2]
    lcut = -17
    ucut = -23
    control = makecontrol(bpt, lier)
    bpt[control] = 8
    maxbpt = int(np.max(bpt))
    nbins = 6
    bins = np.linspace(-23,-17,nbins+1)
    binsize = bins[1] - bins[0]

    #more useful arrays
    #c = make_cmap(int(np.max(bpt) + 1), 'gnuplot')
    c = ['k', 'b', 'r', 'gold', 'g', 'm','g','g','m','b']
    popts = np.zeros((len(Re), maxbpt+1))
    pcovs = np.zeros((len(Re), maxbpt+1))

    plt.figure(figsize=(12,8))
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
            ad = ad/(harc**2)
            ade = np.sqrt((ade/(harc**2))**2 + ((2*ad*harce)/(harc**3))**2)
        
        cut = (Mi < lcut) * (Mi > ucut)

        poptall,pcovall = curve_fit(exponential, Mi[cut], ad[cut], 
                maxfev = 10000, sigma = ade[cut])
        poptall = np.log10(poptall)

        plt.subplot(231+j)
        for k in interested:
            c1 = np.delete((bpt==k)[spxtolier], bad)
            Micut = Mi[c1]
            Miecut = Mie[c1]
            adcut = ad[c1]
            adecut = ade[c1]

            cut = (Mi[c1] < lcut) * (Mi[c1] > ucut)
            adcut = adcut[cut]
            adecut = adecut[cut]
            Micut = Micut[cut]
            Miecut = Miecut[cut]

            means = np.zeros(nbins)
            stds = np.zeros(nbins)
            for l in range(nbins):
                adbin = adcut[(Micut > bins[l]) * (Micut < bins[l+1])]
                adebin = adecut[(Micut > bins[l]) * (Micut < bins[l+1])]
                means[l] = np.average(adbin)
                means[l], stds[l] = weighted_avg_std(adbin, adebin)

            error = stds/means
            plt.errorbar(Micut, np.log10(adcut), yerr = adecut/adcut,
                    fmt='.', c=c[k], elinewidth = .5, ms=3,
                    ecolor = '.75', alpha = .15)
            plt.plot(bins[:-1] + binsize/2, np.log10(means), '-', color=c[k],
                    label=types[k])
            plt.fill_between(bins[:-1] + binsize/2, np.log10(means)+error,
                    np.log10(means) - error, color = c[k], alpha = .2)
            popt,pcov = curve_fit(exponential, Micut, adcut,
                    sigma=adecut, maxfev= 10000)
            popt = np.log10(popt)

            popts[j,k] = popt[1]
            pcovs[j,k] = pcov[1,1]
            plt.legend()
            plt.xlabel(r'$M_i$')
            plt.ylabel(r'$\log AD^2$')
            ax = plt.gca()
            ax.set_ylim((1,6))
            if adr:
                plt.ylabel(r'$\log (AD^2/H_{rot}^2)$')
                ax.set_ylim((-2,0))
            plt.title(r'%s $R_e$' % Re[j])
            plt.grid(True)
        plt.tight_layout()

    '''
    for q in interested:
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
        ax.set_ylim((-1,.25))
    '''
    plt.tight_layout()
    plt.show()
