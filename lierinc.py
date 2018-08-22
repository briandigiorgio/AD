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
    adr = False #set to normalize by gas rotation speed

    #make useful arrays
    spx = np.asarray([spx25, spx50, spx75, spx10, spx12])
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified', 'cLIER below cut', 'SF below cut')
    x = np.linspace(-23,-17,100)
    nbins = 6
    bins = np.linspace(8.5,11.5,nbins+1)#-23,-17,nbins+1)
    binsize = bins[1] - bins[0]

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    maxbpt = int(np.max(bpt))

    #parameters of interest
    interested = [1,2]
    lcut = 12
    ucut = 8
    median = False
    f = False
    incs = np.asarray([])

    if 6 in interested:
        r = 3
        dip = -.6
        group = 2
        ccut(bpt, dip, group, lier, spx, r, i=6)

    if 7 in interested:
        r = 3
        dip = -.6
        group = 1
        ccut(bpt, dip, group, lier, spx, r, i=7)

    #more useful arrays
    #c = make_cmap(int(np.max(bpt) + 1), 'gnuplot')
    c = ['k', 'b', 'r', 'gold', 'g', 'm', 'm', 'g']
    popts = np.zeros((len(Re), maxbpt+1))
    pcovs = np.zeros((len(Re), maxbpt+1))

    for j in [3]:
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
        inc = spx[j]['sini']#['inc_phot']
        ince = spx[j]['sini_e']
        harc = spx[j]['gas_vrot']/inc#['harc_em']
        #harce = spx[j]['harc_se']
        harce = np.sqrt((spx[j]['gas_vrot_e']/inc)**2 +\
                (spx[j]['gas_vrot']*ince/(inc**2))**2)
        #Mi = spx[j]['elpetro_absmag'][:,5]
        #Mie= spx[j]['elpetro_abmerr'][:,5]
        Mi = np.log10(spx[j]['sersic_mass'])
        Mie = np.zeros_like(Mi)
        bad = np.where(np.isnan(np.log(harc*harce*ad*ade)))

        ad = np.delete(ad, bad)
        ade = np.delete(ade, bad)
        harc = np.delete(harc, bad)
        harce = np.delete(harce, bad)
        Mi = np.delete(Mi, bad)
        Mie = np.delete(Mie, bad)
        inc = np.delete(inc, bad)

        for k in interested:
            if not f:
                f, ax1 = plt.subplots(nrows = 2, ncols = 1, figsize = (5,6))
                ax1 = plt.subplot2grid((5,1), (0,0), rowspan = 4)
                ax2 = plt.subplot2grid((5,1), (4,0))

            c1 = np.delete((bpt==k)[spxtolier], bad)
            Micut = Mi[c1]
            Miecut = Mie[c1]
            harccut = harc[c1]
            harcecut = harce[c1]
            inccut = inc[c1]

            cut = (Mi[c1] < lcut) * (Mi[c1] > ucut)# * (inc[c1] > .707)
            harccut = harccut[cut]
            harcecut = harcecut[cut]
            Micut = Micut[cut]
            Miecut = Miecut[cut]
            inccut = inccut[cut]

            means = np.zeros(nbins)
            stds = np.zeros(nbins)
            for l in range(nbins):
                harcbin = harccut[(Micut > bins[l]) * (Micut < bins[l+1])]
                harcebin = harcecut[(Micut > bins[l]) * (Micut < bins[l+1])]
                means[l], stds[l] = weighted_avg_std(harcbin, harcebin)

                if median:
                    means[l] = np.median(harcbin)

            error = stds
            #plt.errorbar(Micut, harccut, yerr = harcecut,
            #        fmt='.', elinewidth = .5, ms=3,
            #        ecolor = '.75')
            incp = ax1.scatter(Micut, harccut, c=inccut, cmap = 'gnuplot', s=20)
            #ax1.plot(bins[:-1] + binsize/2, means, '-', color=c[k],
            #        label=types[k])
            #ax1.fill_between(bins[:-1] + binsize/2, means + error,
            #        means - error, color = c[k], alpha = .2)

            ax1.legend()
            #cbar = plt.colorbar(inc)
            #cbar.set_label('Inclination (degrees)')
            ax1.set_xlabel(r'$M_i$')
            ax1.set_ylabel(r'$V_c$')
            ax1.set_ylim((40,400))
            ax1.set_xlim((8.5,11.5))#-17,-23))
            ax1.set_title(r'Tully-Fisher at %s $R_e$' % Re[j])
            ax1.grid(True)
            ax1.set_yscale('log')

            ax2.set_xlabel(r'$\sin$(Inclination)')
            ax2.set_xlim((0,1))#90))
            incs=np.append(incs,inccut)

        colorhist(incs, 'gnuplot', bins = 20, ax = ax2)
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
