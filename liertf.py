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
    summary = np.load('modelfitting/summary.npy')
    spx = readall()
    drp = loadfits('drpall-v2_3_1.fits')
    mass = True
    masstype = 'nsa_elpetro_mass' #drpall mass variable
    conc = 'pet' #pet = pet90/pet50, ser = sersic n

    drpplate = drp['plate']
    drpifu = [int(i) for i in drp['ifudsgn']]
    lierplate = lier['PLATE']
    lierifu = [int(i) for i in lier['IFUDESIGN']]

    liertodrp = plateifu(drpplate, drpifu, lierplate, lierifu)

    #make useful arrays
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified', 'cLIER below cut', 'SF below cut', 'SF Control')
    c = ['k', 'b', 'r', 'gold', 'g', 'm', 'm', 'g','b']

    #parameters of interest
    interested = [2,8]
    lcut = -17
    ucut = -23
    if mass:
        lcut = 8.5
        ucut = 11.5
    median = False
    nbins = 6
    bins = np.linspace(ucut,lcut,nbins+1)
    if mass:
        bins = np.linspace(lcut,ucut,nbins+1)
    binsize = bins[1] - bins[0]

    #for j in [3]:
    j=3

    #get plate/ifu data for matching
    plate = spx[j]['plate']
    ifu = spx[j]['ifudesign']
    plateifuspx = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])
    spxtodrp = plateifu(drpplate,drpifu,plate,ifu)
    spxtolier = plateifu(lierplate,lierifu,plate,ifu)

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    control = makecontrol(bpt, lier, mass = drp[masstype][liertodrp],
            c=conc)#, stl = spxtolier, spx = spx[j])
    #control = makecontrol(bpt, lier, c=conc)
    bpt[control] = 8
    maxbpt = int(np.max(bpt))

    #filter out bad values and pick out correct data
    ad = spx[j]['ad2_em']
    ade =spx[j]['ad2_se']
    inc = spx[j]['sini']
    ince = spx[j]['sini_e']
    #harc = spx[j]['harc_em']
    harc = spx[j]['gas_vrot']/inc
    #harce = spx[j]['harc_se']
    harce = np.sqrt((spx[j]['gas_vrot_e']/inc)**2 +\
            (spx[j]['gas_vrot']*ince/(inc**2))**2)
    Mi = spx[j]['elpetro_absmag'][:,5]
    Mie= spx[j]['elpetro_abmerr'][:,5]
    if mass:
        Mi = np.log10(drp[masstype][spxtodrp])
        Mie = np.zeros_like(Mi)
    bad = np.where(np.isnan(np.log(harc*harce*ad*ade)))

    spxtolier = np.delete(spxtolier, bad)
    ad = np.delete(ad, bad)
    ade = np.delete(ade, bad)
    harc = np.delete(harc, bad)
    harce = np.delete(harce, bad)
    Mi = np.delete(Mi, bad)
    Mie = np.delete(Mie, bad)

    for k in range(len(interested)):
    #for k in interested:
        c1 = (bpt==interested[k])[spxtolier]
        Micut = Mi[c1]
        Miecut = Mie[c1]
        harccut = harc[c1]
        harcecut = harce[c1]

        #plt.hist(np.delete(spx[j]['reff'], bad)[c1], histtype = 'step', 
        #        label = types[interested[k]], bins = np.linspace(0,20,20))

        cut = (Mi[c1] < lcut) * (Mi[c1] > ucut) * (inc[c1] > .707)
        if mass:
            cut = (Mi[c1] > lcut) * (Mi[c1] < ucut) * (inc[c1] > .707)
        harccut = harccut[cut]
        harcecut = harcecut[cut]
        Micut = Micut[cut]
        Miecut = Miecut[cut]

        means = np.zeros(nbins)
        stds = np.zeros(nbins)
        for l in range(nbins):
            harcbin = harccut[(Micut > bins[l]) * (Micut < bins[l+1])]
            harcebin = harcecut[(Micut > bins[l]) * (Micut < bins[l+1])]
            means[l], stds[l] = weighted_avg_std(harcbin, harcebin)

            if median:
                means[l] = np.median(harcbin)

            '''
            weird = (np.abs(harcbin-means[l]) > 75)
            wp = plateifuspx[c1][cut][(Micut > bins[l]) * \
                    (Micut < bins[l+1])][weird]
            wv = harcbin[weird] - means[l]
            print(types[interested[k]], bins[l], '-', bins[l+1])
            print(*list(sorted(zip(wp,wv), key=lambda t: -t[1])), sep='\n')
            print()
            '''

        error = stds
        print('Number of galaxies of type %d: %d' % (interested[k],len(Micut)))
        plt.errorbar(Micut, harccut, yerr = harcecut,
                fmt='.', c=c[interested[k]], elinewidth = .5, ms=5,
                ecolor = '.75')
        plt.plot(bins[:-1] + binsize/2, means, '-', color=c[interested[k]],
                label=types[interested[k]])
        plt.fill_between(bins[:-1] + binsize/2, means + error,
                means - error, color = c[interested[k]], alpha = .2)

        plt.legend()
        plt.xlabel(r'$M_i$')
        if mass:
            plt.xlabel('log ' + masstype)
        plt.ylabel(r'$V_{rot}$')
        ax = plt.gca()
        ax.set_ylim((0,400))
        ax.set_xlim((lcut,ucut))
        #ax.set_yscale('log')
        plt.title(r'Tully-Fisher, concentration = ' + conc)
        plt.grid(True)
    plt.tight_layout()
    
    for k in range(len(types)):
        kcut = (bpt==k)[spxtolier]
        plt.plot(Mi[kcut], harc[kcut], '.', c = c[k], alpha = .2, ms = 3)


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
