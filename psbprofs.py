#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *
import time
import matplotlib.lines as mlines

if __name__ == '__main__':
    #load all relevant files
    start = time.time()
    match = loadfits('manga_catalog_match.fits.gz', i=2)
    spx = readall()
    adr = False #set to normalize by gas rotation speed
    drpall = loadfits('drpall-v2_3_1.fits')

    #make useful arrays
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    c = ['k', 'b', 'r', 'gold', 'g', 'm', 'b']

    #parameters of interest
    interested = [0,1]
    lcut = 11.5
    ucut = 8.5
    median = False
    nbins = 3
    bins = np.linspace(ucut,lcut,nbins+1)
    binsize = bins[1] - bins[0]

    #more useful arrays
    #c = make_cmap(int(np.max(bpt) + 1), 'gnuplot')
    adprof  = np.zeros((len(Re), len(interested), nbins))
    adeprof = np.zeros((len(Re), len(interested), nbins))
    harcprof  = np.zeros((len(Re), len(interested), nbins))
    harceprof = np.zeros((len(Re), len(interested), nbins))
    strcprof  = np.zeros((len(Re), len(interested), nbins))
    strceprof = np.zeros((len(Re), len(interested), nbins))

    for j in range(len(Re)):
        #get plate/ifu data for matching
        plate = spx[j]['plate']
        ifu = spx[j]['ifudesign']
        #plateifuspx = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])

        drpplate = drpall['plate']
        drpifu = np.asarray([int(i) for i in drpall['ifudsgn']])
        #drpifu = drpall['ifudsgn'].astype(str)
        #plateifudrp = np.asarray([drpplate[i] + drpifu[i] 
        #    for i in range(len(drpplate))])

        #match the catalogs
        #for some reason the numpy version doesn't work but the python one does
        #spxtolier = np.asarray([np.argmax(plateifulier==plateifuspx[i])
        #                for i in range(len(plateifuspx))])

        #spxtodrp = np.zeros(len(plateifuspx))
        #for l in range(len(plateifuspx)):
        #    for m in range(len(plateifudrp)):
        #        if plateifuspx[l] == plateifudrp[m]:
        #            spxtodrp[l] = m
        #spxtodrp = spxtodrp.astype(int)
        spxtodrp = plateifu(drpplate, drpifu, plate, ifu)

        #filter out bad values and pick out correct data
        ad  = spx[j]['ad2_em']
        ade = spx[j]['ad2_se']
        harc  = spx[j]['harc_em']
        harce = spx[j]['harc_se']
        strc  = spx[j]['strc_em']
        strce = spx[j]['strc_se']
        #Mi  = spx[j]['elpetro_absmag'][:,5]
        #Mie = spx[j]['elpetro_abmerr'][:,5]
        Mi = np.log10(spx[j]['sersic_mass'])
        Mie = np.zeros_like(Mi)
        bad = np.where(np.isnan(np.log(harc*harce*ad*ade*strc*strce)))

        spxtodrp = np.delete(spxtodrp, bad)
        ad  = np.delete(ad, bad)
        ade = np.delete(ade, bad)
        harc  = np.delete(harc, bad)
        harce = np.delete(harce, bad)
        strc  = np.delete(strc, bad)
        strce = np.delete(strce, bad)
        Mi  = np.delete(Mi, bad)
        Mie = np.delete(Mie, bad)

        bulgetypes = np.ones(len(ad))
        bulgetypes[drpall[spxtodrp]['nsa_sersic_n'] > 2] = 0

        #harc[np.isnan(harc) or not harc] = 1
        if adr:
            ad = ad/(harc**2)
            ade = np.sqrt((ade/(harc**2))**2 + ((2*ad*harce)/(harc**3))**2)
        
        cut = (Mi < lcut) * (Mi > ucut)

        for k in range(len(interested)):
            c1 = (bulgetypes==interested[k])
            Micut  = Mi[c1]
            Miecut = Mie[c1]
            adcut  = ad[c1]
            adecut = ade[c1]
            harccut  = harc[c1]
            harcecut = harce[c1]
            strccut  = strc[c1]
            strcecut = strce[c1]
            '''
            adprof[j,k], adeprof[j,k] = weighted_avg_std(adcut, 
                    adecut)
            harcprof[j,k], harceprof[j,k] = weighted_avg_std(harccut, 
                    harcecut)
            strcprof[j,k], strceprof[j,k] = weighted_avg_std(strccut, 
                    strcecut)
            label = ''
            '''
            for l in range(nbins):
                cut = (Mi[c1] < bins[l+1]) * (Mi[c1] > bins[l])
                adcutl  = adcut[cut]
                adecutl = adecut[cut]
                Micutl  = Micut[cut]
                Miecutl = Miecut[cut]
                harccutl  = harccut[cut]
                harcecutl = harcecut[cut]
                strccutl  = strccut[cut]
                strcecutl = strcecut[cut]
        
                adprof[j,k,l], adeprof[j,k,l] = weighted_avg_std(adcutl, 
                        adecutl)
                harcprof[j,k,l], harceprof[j,k,l] = weighted_avg_std(harccutl, 
                        harcecutl)
                strcprof[j,k,l], strceprof[j,k,l] = weighted_avg_std(strccutl, 
                        strcecutl)

                label = ''
                if median:
                    adprof[j,k,l] = np.median(adcutl)
                    harcprof[j,k,l] = np.median(harccutl)
                    strcprof[j,k,l] = np.median(strccutl)
                    label = 'Median '
            
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
    '''
    plt.figure(figsize=(12,4))
    c2 = make_cmap(nbins, 'gnuplot2')
    handles = []
    ls = []

    if 0 in interested:
        solid = mlines.Line2D([],[], color = 'k', label = 'Pseudo')
        handles += [solid]
        ls += ['-']
    if 1 in interested:
        dashed = mlines.Line2D([],[], color = 'k', ls = '--',
            label = 'Classical')
        handles += [dashed]
        ls += ['--']
    if 6 in interested:
        dotted = mlines.Line2D([],[], color = 'k', ls = ':', 
            label = '%s < %g' % (types[group], dip))
        handles += [dotted]
        ls += [':']
    if 8 in interested:
        solid = mlines.Line2D([],[], color = 'k', label = 'SF Control')
        handles += [solid]
        ls += ['-']
    
    print('plotting')
    print(harcprof, strcprof, adprof)
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
                plt.ylabel('%sRotation Speed' % label)
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
