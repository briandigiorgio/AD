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
    lier = loadfits('MPL-6_master_catalogue_Apr9.fits')
    adcut = np.load('modelfitting/adcut.npy')
    hyb = loadfits('adsample.fits')[adcut]
    drp = loadfits('drpall-v2_4_3.fits')

    #parameters
    interested = [2,8]
    masstype = 'nsa_elpetro_mass' #drpall mass variable
    conc = 'ser' #pet = pet90/pet50, ser = sersic n
    qbins = True
    nbins = 6
    median = True
    lcut = 8.5
    ucut = 11.5

    drpplate = drp['plate']
    drpifu = [int(i) for i in drp['ifudsgn']]
    lierplate = lier['PLATE']
    lierifu = [int(i) for i in lier['IFUDESIGN']]
    #liertodrp = plateifu(drpplate, drpifu, lierplate, lierifu)
    #np.save('liertodrp', liertodrp)
    liertodrp = np.load('liertodrp.npy')

    #make useful arrays
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified', 'cLIER below cut', 'SF below cut', 'SF Control')
    c = ['k', 'b', 'r', 'gold', 'g', 'm', 'm', 'g','b']

    #parameters of interest and bin specification
    bins = np.linspace(lcut,ucut,nbins+1)
    binsize = bins[1] - bins[0]

    #get plate/ifu data for matching and match catalogs
    plate = hyb['plate'].astype(int)
    ifu = hyb['ifu'].astype(int)
    #plateifuhyb = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])
    #hybtodrp = plateifu(drpplate,drpifu,plate,ifu)
    #hybtolier = plateifu(lierplate,lierifu,plate,ifu)
    #np.save('hybtodrplier', np.array(list(zip(hybtodrp, hybtolier))).T)
    hybtodrp, hybtolier = np.load('hybtodrplier.npy')

    #clean up the nans in the bpt designations and make SF control
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    control, weights = makecontrol(bpt, lier, mass = drp[masstype][liertodrp],
            c=conc, norm = True, counts = True)
    #np.save('controloutw', np.array(list(zip(control,weights))).T)
    #control, weights = np.load('controloutw.npy')
    bpt[control] = 8

    #put the weights in the right shape to deal with
    wa = np.zeros_like(bpt)
    wa[control] = weights
    wahyb = wa[hybtolier]

    '''
    #look at images for some example galaxies
    img = 2
    imgplate = lier[bpt==img]['PLATE'][:60]
    imgifu = [int(i) for i in lier[bpt==img]['IFUDESIGN'][:60]]
    sdsscollage(imgplate, imgifu, nrow = 10, ncol = 6, npix = 500, scale = .2,
            drpall = drp)
    plt.show()
    '''

    #load relevant parameters and errors
    #inc = np.sin(hyb['inc'] * np.pi/180)
    #ince = np.cos(hyb['inc'] * np.pi/180) * hyb['inc_err'] * np.pi/180
    #harc = hyb['gas_vrot']/inc
    #harce = np.sqrt((hyb['gas_vrot_err']/inc)**2 +\
    #        (hyb['gas_vrot']*ince/(inc**2))**2)
    inc = np.sin(hyb['inc'] * np.pi/180)
    harc = hyb['vmod']/inc
    harce = hyb['vmode']/inc
    #harc = hyb['svrot'][:,14]/inc
    #harce = hyb['svrote'][:,14]/inc

    Mi = np.log10(drp[masstype][hybtodrp])
    Mie = np.zeros_like(Mi)

    #make correct number of bins with equal numbers of galaxies
    #if qbins:
    #    bincut = np.zeros_like(Mi)
    #    for k in range(len(interested)):
    #        bincut += (bpt==interested[k])[hybtolier].astype(bool)
    #    bincut = bincut.astype(bool)
    #    bins = np.percentile(Mi[bincut], np.linspace(0,100,nbins+1))
    bins = makeqbins(nbins, Mi, bpt, hybtolier, interested)

    print('Bins: ',bins)

    #for all interested bpt types, plot TF and average in each bin
    for k in range(len(interested)):
        #cut data to only be interested bpt type
        c1 = (bpt==interested[k])[hybtolier]
        Micut = Mi[c1]
        Miecut = Mie[c1]
        harccut = harc[c1]
        harcecut = harce[c1]
        weightscut = wahyb[c1]

        #cut based off magnitude and inclination
        cut = (Mi[c1] > lcut) * (Mi[c1] < ucut) * (inc[c1] > .707)
        harccut = harccut[cut]
        harcecut = harcecut[cut]
        Micut = Micut[cut]
        Miecut = Miecut[cut]
        weightscutl = weightscut[cut]
        print('Number of galaxies of type %d: %d' % (interested[k],len(Micut)))

        #calculate means and stds for each bin of data
        means = np.zeros(nbins)
        stds = np.zeros(nbins)
        for l in range(nbins):
            harcbin = harccut[(Micut > bins[l]) * (Micut < bins[l+1])]
            harcebin = harcecut[(Micut > bins[l]) * (Micut < bins[l+1])]
            if interested[k] == 8:
                weightsbin = weightscutl[(Micut > bins[l])*(Micut < bins[l+1])]
                means[l], stds[l] = weighted_avg_std(harcbin, 1/weightsbin)
            else:
                #harcebin = np.ones_like(harcebin)
                means[l], stds[l] = weighted_avg_std(harcbin, harcebin)

            if median:
                means[l] = np.median(harcbin)

        plt.errorbar(Micut, harccut, yerr = harcecut, fmt='.', ms = 5,
                c = c[interested[k]], elinewidth = .5, ecolor = '.75')
        plt.semilogy(bins[:-1] + .5*(bins[1:]-bins[:-1]), means, '-o', 
                color = c[interested[k]], label = types[interested[k]])
        plt.fill_between(bins[:-1] + .5 * (bins[1:] - bins[:-1]), means + stds,
                means - stds, color = c[interested[k]], alpha = .2)

        plt.legend()
        plt.xlabel('log ' + masstype)
        plt.ylabel(r'$V_{rot}$')
        ax = plt.gca()
        ax.set_ylim((0,400))
        ax.set_xlim((lcut,ucut))
        plt.title(r'Tully-Fisher, concentration = ' + conc)
        #plt.grid(True)

    plt.tight_layout()
    [plt.axvline(b, alpha = .15, color = 'k', ls = '--') for b in bins]
    
    #plot the rest of the unused galaxies in the ad sample in the background
    for k in range(len(types)):
        kcut = (bpt==k)[hybtolier]
        plt.plot(Mi[kcut], harc[kcut], '.', c = c[k], alpha = .15, ms = 3)

    plt.tight_layout()
    plt.show()
