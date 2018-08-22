#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *
import time
import matplotlib.lines as mlines

if __name__ == '__main__':
    #load all relevant files
    start = time.time()
    lier = loadfits('gal_list_v2_0_1_bpt_classify3.fits')
    match = loadfits('manga_catalog_match.fits.gz', i=2)
    spx = readall()
    adr = False #set to normalize by gas rotation speed

    #make useful arrays
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified', 'SF below cut', 'cLIER below cut')
    x = np.linspace(-23,-17,100)
    nbins = 3
    bins = np.linspace(-23,-17,nbins+1)
    binsize = bins[1] - bins[0]

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    maxbpt = int(np.max(bpt))

    #parameters of interest
    interested = [1, 2, 6,7]
    lcut = -21
    ucut = -23

    #cut the cLIERs into two populations as seen in lierhist
    if 6 in interested:
        r = 3 #radius in Re list to make cut on 
        dip = -.6 #dip in bimodal histogram to cut on
        ccutul(bpt, dip, 1, lier, spx, r, i = 6, l = -21, u = -23)
        ccutul(bpt, dip, 2, lier, spx, r, i = 7, l = -21, u = -23)

    #more useful arrays
    #c = make_cmap(int(np.max(bpt) + 1), 'gnuplot')
    c = ['k', 'b', 'r', 'gold', 'g', 'm']
    cmaps = ['Blues', 'Reds', 'RdPu', 'Greens']

    adcut  = [[]] * len(Re)
    adecut = [[]] * len(Re)
    harccut  = [[]] * len(Re)
    harcecut = [[]] * len(Re)
    strccut  = [[]] * len(Re)
    strcecut = [[]] * len(Re)
    '''
    adcut  = np.zeros((len(Re), 607))
    adecut = np.zeros((len(Re), 607))
    harccut  = np.zeros((len(Re), 607))
    harcecut = np.zeros((len(Re), 607))
    strccut  = np.zeros((len(Re), 607))
    strcecut = np.zeros((len(Re), 607))
    '''
    spxtolier = [[]] * len(Re)

    #get plate/ifu data for matching
    j=3

    lierplate = lier['PLATE'].astype(str)
    lierifu = lier['IFUDESIGN'].astype(str)
    plateifulier = np.asarray([lierplate[i] + lierifu[i] 
        for i in range(len(lierplate))])

    #match the catalogs
    #for some reason the numpy version doesn't work but the python one does

    for b in range(len(Re)):
        plate = spx[b]['plate'].astype(str)
        ifu = spx[b]['ifudesign'].astype(str)
        plateifuspx = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])
        spxtolier[b] = np.zeros(len(plateifuspx))
        #spxtolier = np.asarray([np.argwhere(plateifulier==plateifuspx[i])
        #            for i in range(len(plateifuspx))])
        for a in range(len(plateifuspx)):
            for m in range(len(plateifulier)):
                if plateifuspx[a] == plateifulier[m]:
                    spxtolier[b][a] = m
        spxtolier[b] = np.asarray(spxtolier[b]).astype(int)
    spxtolier = np.asarray(spxtolier)

    #filter out bad values and pick out correct data
    ad  = spx[j]['ad2_em']
    ade = spx[j]['ad2_se']
    harc  = spx[j]['harc_em']
    harce = spx[j]['harc_se']
    strc  = spx[j]['strc_em']
    strce = spx[j]['strc_se']
    Mi  = spx[j]['elpetro_absmag'][:,5]
    Mie = spx[j]['elpetro_abmerr'][:,5]
    bad = np.where(np.isnan(np.log(harc*harce*ad*ade*strc*strce)))
    '''
    spxtolier = np.delete(spxtolier, bad)
    ad  = np.delete(ad, bad)
    ade = np.delete(ade, bad)
    harc  = np.delete(harc, bad)
    harce = np.delete(harce, bad)
    strc  = np.delete(strc, bad)
    strce = np.delete(strce, bad)
    Mi  = np.delete(Mi, bad)
    Mie = np.delete(Mie, bad)
    '''
    #harc[np.isnan(harc) or not harc] = 1
    if adr:
        ad = ad/(harc**2)
        ade = np.sqrt((ade/(harc**2))**2 + ((2*ad*harce)/(harc**3))**2)
    
    plotnum = 0
    plt.figure(figsize=(8,8))
    for k in range(len(interested)):
        for l in range(len(Re)):
            cut = (bpt==interested[k])[spxtolier[l]]
            print(interested[k], Re[l], np.sum(cut))
            adcut[l]  = spx[l]['ad2_em'][cut]
            adecut[l] = spx[l]['ad2_se'][cut]
            harccut[l]  = spx[l]['harc_em'][cut]
            harcecut[l] = spx[l]['harc_se'][cut]
            strccut[l]  = spx[l]['strc_em'][cut]
            strcecut[l] = spx[l]['strc_se'][cut]
        
        cmap = make_cmap(len(adcut[l]), cmaps[k])

        plt.subplot(221+plotnum)
        for m in range(np.min([len(adcut[d]) for d in range(len(Re))])):
            data = [harccut[d][m] for d in range(len(Re))]
            plt.plot(Re, data, c = cmap[m])
        plt.grid(True)
        plt.title(types[interested[k]])
        plt.xlabel(r'$R_e$')
        plt.ylabel(r'$V_c$')
        ax = plt.gca()
        ax.set_ylim((0,400))
        plotnum += 1
    plt.tight_layout()
    plt.show()

    '''
    plt.figure(figsize=(12,8))
    for l in range(nbins):
        plt.subplot(231+l)
        for k in range(len(interested)):
            plt.plot(Re, np.sqrt(adprof[:,k,l]), '-', color = c[interested[k]],
                    label = 'AD for %s' % types[interested[k]])
            plt.plot(Re, harcprof[:,k,l], '--', color = c[interested[k]],
                    label = 'Gas for %s' % types[interested[k]])
            plt.plot(Re, strcprof[:,k,l], ':', color = c[interested[k]],
                    label = 'Stars for %s' % types[interested[k]])
            

        plt.title(r'$M_i = $(%d,%d)' % (bins[l+1], bins[l]))
        plt.legend()
        plt.grid(True)

    plt.figure(figsize=(12,4))
    c2 = make_cmap(nbins, 'gnuplot2')
    handles = []
    ls = []

    if 1 in interested:
        solid = mlines.Line2D([],[], color = 'k', label = 'Star-Forming')
        handles += [solid]
        ls += ['-']
    if 2 in interested:
        dashed = mlines.Line2D([],[], color = 'k', ls = '--',
            label = 'cLIER')
        handles += [dashed]
        ls += ['--']
    if 6 in interested:
        dotted = mlines.Line2D([],[], color = 'k', ls = ':', 
            label = '%s < %g' % (types[group], dip))
        handles += [dotted]
        ls += [':']

    for k in range(len(interested)):
        for m in range(3):
            for l in range(nbins):
                plt.subplot(131)
                plt.plot(Re, harcprof[:,k,l], ls[k], color = c2[-l-1])
                plt.fill_between(Re, harcprof[:,k,l]-harceprof[:,k,l],
                        harcprof[:,k,l]+harceprof[:,k,l], color = c2[-l-1],
                        alpha = .07)
                plt.title('Gas Rotation Speed')
                plt.xlabel(r'$R_e$')
                plt.ylabel('Rotation Speed')
                ax = plt.gca()
                ax.set_ylim((0, 400))
                plt.grid(True)

                plt.subplot(132)
                plt.plot(Re, strcprof[:,k,l], ls[k], color = c2[-l-1])
                plt.fill_between(Re, strcprof[:,k,l]-strceprof[:,k,l],
                        strcprof[:,k,l]+strceprof[:,k,l], color = c2[-l-1],
                        alpha = .07)
                plt.title('Stellar Rotation Speed')
                plt.xlabel(r'$R_e$')
                ax = plt.gca()
                ax.set_ylim((0, 400))
                ax.tick_params(left = False, labelleft = False)
                plt.grid(True)

                plt.subplot(133)
                plt.plot(Re, np.sqrt(adprof[:,k,l]), ls[k], color = c2[-l-1])
                plt.fill_between(Re, np.sqrt(adprof[:,k,l]-adeprof[:,k,l]),
                        np.sqrt(adprof[:,k,l]+adeprof[:,k,l]), color = c2[-l-1],
                        alpha = .07)
                plt.title('AD')
                plt.xlabel(r'$R_e$')
                ax = plt.gca()
                ax.set_ylim((0, 400))
                ax.tick_params(left = False, labelleft = False)
                plt.grid(True)



    plt.legend(handles = handles)
    plt.tight_layout()
    plt.show()
    '''
