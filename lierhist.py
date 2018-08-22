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
    mass = True

    #make useful arrays
    spx = np.asarray([spx25, spx50, spx75, spx10, spx12])
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified', 'cLIER cut','','SF control')
    
    x = np.linspace(-23,-17,100)

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    maxbpt = int(np.max(bpt))
    control = makecontrol(bpt, lier)
    bpt[control] = 8
    maxbpt = int(np.max(bpt))

    #parameters of interest
    interested = [1,2]
    lcut = -17
    ucut = -23
    if mass:
        lcut = 11.5
        ucut = 8.5
    nbins = 6
    bins = np.linspace(ucut,lcut,nbins+1)
    binsize = bins[1] - bins[0]

    #more useful arrays
    #c = make_cmap(int(np.max(bpt) + 1), 'gnuplot')
    c = ['k', 'b', 'r', 'gold', 'g', 'm', 'brown','g','b']
    popts = np.zeros((len(Re), maxbpt+1))
    pcovs = np.zeros((len(Re), maxbpt+1))

    #cut the cLIERs into two populations
    if 6 in interested:
        r = 3 #radius in Re list to make cut on 
        dip = -.6 #dip in bimodal histogram to cut on
        group = 2 #type number of galaxies to cut (1=sf, 2 = clier)
        plate = spx[r]['plate'].astype(str)
        ifu = spx[r]['ifudesign'].astype(str)
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

        ad10 = spx[r]['ad2_em']
        harc10  = spx[r]['harc_em']
        bpt[spxtolier * (bpt[spxtolier]==group) * \
                (np.log10(ad10/(harc10**2)) < dip)] = 6

    for j in range(len(Re)):
        f, ax1 = plt.subplots(2, ncols = 2, sharey = True)
        ax1 = plt.subplot2grid((1,4), (0,0), colspan = 3)
        ax2 = plt.subplot2grid((1,4), (0,3))
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
        if mass:
            Mi = np.log10(spx[j]['sersic_mass'])
            Mie = np.zeros_like(Mi)
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

            ax1.errorbar(Micut, np.log10(adcut), yerr = adecut/adcut,
                    fmt='.', c=c[k], elinewidth = .5, ms=3,
                    ecolor = '.75', alpha = .25)
            ax1.plot(bins[:-1] + binsize/2, np.log10(means), '-', color=c[k],
                    label=types[k])
            ax1.fill_between(bins[:-1] + binsize/2, np.log10(means)+error,
                    np.log10(means) - error, color = c[k], alpha = .2)
            popt,pcov = curve_fit(exponential, Micut, adcut,
                    sigma=adecut, maxfev= 10000)
            popt = np.log10(popt)

            popts[j,k] = popt[1]
            pcovs[j,k] = pcov[1,1]
            ax1.legend()
            ax1.set_xlabel(r'$M_i$')
            if mass:
                ax1.set_xlabel('log NSA Sersic Mass')
            ax1.set_ylabel(r'$\log AD^2$')
            ax1.set_ylim((0,6))
            ax2.set_ylim((0,6))
            ax1.grid(True)
            if adr:
                ax1.set_ylabel(r'$\log (AD^2/V_c^2)$')
                ax1.set_ylim((-2,0))
                ax2.set_ylim((-2,0))
            ax1.set_title(r'%s $R_e$' % Re[j])
            ax1.grid(True)

            ax2.hist(np.log10(adcut), orientation = 'horizontal', color = c[k],
                    density = True, alpha = .4, bins = 25, range = (-2,0))
            ax2.grid(True)
            ax2.tick_params('both', left = False, labelleft = False, bottom
                    = False, labelbottom = False)
            f.subplots_adjust(wspace=0.05)
    plt.show()
