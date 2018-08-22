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
    adr = False #do AD/gas rotation speed^2 instead

    #make useful arrays
    spx = np.asarray([spx25, spx50, spx75, spx10, spx12])
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified')
    x = np.linspace(-23,-17,100)

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    maxbpt = int(np.max(bpt))
    lcut = -17
    ucut = -23
    interested = [1] #[0,1,2,3,5]

    #more useful arrays
    #c = make_cmap(int(np.max(bpt) + 1), 'gnuplot')
    c = ['k', 'b', 'r', 'gold', 'g', 'm']
    popts = np.zeros(len(Re))
    pcovs = np.zeros(len(Re))
    poptsall = np.zeros(len(Re))
    pcovsall = np.zeros(len(Re))

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

        #get all relevant data
        ad = spx[j]['ad2_em']
        ade =spx[j]['ad2_se']
        harc = spx[j]['harc_em']
        harce = spx[j]['harc_se']
        Mi = spx[j]['elpetro_absmag'][:,5]
        Mie= spx[j]['elpetro_abmerr'][:,5]

        #pick out bad values (negatives, nans, 0s, etc.)
        bad = np.where(np.isnan(np.log(harc*harce*ad*ade)))
        ad = np.delete(ad, bad)
        ade = np.delete(ade, bad)
        harc = np.delete(harc, bad)
        harce = np.delete(harce, bad)
        Mi = np.delete(Mi, bad)
        Mie = np.delete(Mie, bad)

        #changes to ad/rotation speed if specified
        if adr:
            ad = ad/(harc**2)
            ade = np.sqrt((ade/(harc**2))**2 + ((2*ad*harce)/(harc**3))**2)
        ''' 
        #make plot for comparing all of the different types
        plt.subplot(326)
        plt.title('Fit for All BPT Classifications')
        plt.xlabel(r'$M_i$')
        plt.ylabel(r'$\Delta AD^2$')
        plt.grid(True)
        ax = plt.gca()
        ax.set_ylim((-1.5,1.5))
        plt.axhline(y=0, color = 'k')
        '''
        #cut for looking at only a part of the magnitudes
        cut = (Mi < lcut) * (Mi > ucut)

        #fitline for all of the data
        poptall,pcovall = curve_fit(exponential, Mi[cut], ad[cut],
                maxfev=10000, sigma=ade[cut])
        poptall = np.log10(poptall) #transform to linear
        poptsall[j] = poptall[1]
        pcovsall[j] = pcovall[1,1]

        plt.subplot(321+j)
        for k in interested:
            #match and cut bpt data
            c1 = np.delete((bpt==k)[spxtolier], bad)
            Micut = Mi[c1]
            Miecut = Mie[c1]
            adcut = ad[c1]
            adecut = ade[c1]

            #magnitude cut
            cut = (Mi[c1] < lcut) * (Mi[c1] > ucut)
            adcut = adcut[cut]
            adecut = adecut[cut]
            Micut = Micut[cut]
            Miecut = Miecut[cut]

            error = adecut/adcut #log error for plots, not for fit line

            #plot difference between fit and data for each type
            #plt.subplot(321+k)

            #find fit line by fitting in log space and transforming to linear
            popt,pcov = curve_fit(exponential, Micut, adcut,
                    maxfev=100000, sigma=adecut)
            popt = np.log10(popt)
            popts[j] = popt[1]
            pcovs[j] = pcov[1,1]

            #plot difference between fit and data
            delta = np.log10(adcut) - line(Micut, popt[1], popt[0])
            plt.errorbar(Micut, delta, xerr=Miecut,yerr = error, fmt='.', 
                    c=c[k], label = types[k], elinewidth = .5, ms = 3, 
                    ecolor='.75')
            
            #plot overall fitline on graph too
            plt.plot(x, line(x, poptall[1], poptall[0]) - 
                    line(x, popt[1], popt[0]), 'k--')
            plt.axhline(y=0, color = 'k')

            #plt.legend()
            plt.xlabel(r'$M_i$')
            plt.ylabel(r'$\Delta \log AD^2$')
            ax = plt.gca()
            ax.set_xlim((-23,-17))
            ax.set_ylim((-1.5,1.5))
            if adr:
                plt.ylabel(r'$\Delta \log AD^2/H_{rot}^2$')
            #plt.title(types[k])
            plt.title(r'%s $R_e$' %Re[j])
            plt.grid(True)
            plt.axhline(y=0, color = 'k')
            '''
            #plot diff between overall fit and data on separate plot
            plt.subplot(326)
            deltall = np.log10(adcut) - line(Micut, poptall[1], poptall[0])
            plt.errorbar(Micut, deltall, xerr=Miecut, yerr = error,fmt='.', 
                    c=c[k], label = types[k], elinewidth = .5, ms = 3)
            #plot all of the other fit lines, makes plot pretty messy
            #plt.plot(x, line(x,popt[1], popt[0]) - 
            #    line(x, poptall[1], poptall[0]), color = c[k])
            '''
        plt.tight_layout()

    plt.subplot(326)
    plt.title('Slopes for All Radii')
    plt.xlabel(r'$R_e$')
    plt.ylabel('Slope')
    plt.grid(True)
    ax = plt.gca()
    #ax.set_ylim((-1.5,1.5))
    #plt.axhline(y=0, color = 'k')
    print(popts,pcovs)
    plt.plot(Re, popts, color='k', label = types[interested[0]])
    plt.fill_between(Re, popts-pcovs, popts+pcovs, color='k',
            alpha = .2)
    plt.plot(Re, poptsall, 'k--', label = 'All')
    plt.legend()
    plt.show()
