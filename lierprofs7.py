#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *
import time
import matplotlib.lines as mlines
from astropy import constants as const

if __name__ == '__main__':
    #load all relevant files
    start = time.time()
    lier = loadfits('MPL-6_master_catalogue_Apr9.fits')
    drp = loadfits('drpall-v2_4_3.fits')
    hyb = loadfits('adsample.fits')
    bigcut = hyb['cut']
    hyb = hyb[bigcut]
    conc = 'ser'
    masstype = 'MASS_SER'#'MASS_ELL_PETRO'
    hybtodrp, hybtolier = np.load('hybtodrplier.npy')
    liertodrp = np.load('liertodrp.npy')

    ad2 = hyb['ad']
    ad2e = hyb['ade']
    gvrot = hyb['gvrot']
    gvrote = hyb['gvrote']
    svrot = hyb['svrot']
    svrote = hyb['svrote']
    plate = hyb['plate']
    ifu = hyb['ifu']
    Mi = lier[masstype][hybtolier]
    Mie = np.zeros_like(Mi)
    re = np.load('re.npy')[bigcut][:,0]

    #match lier with data
    lierplate = lier['PLATE']
    lierifu = np.asarray([int(i) for i in lier['IFUDESIGN']])
    hybtolier = plateifu(lierplate,lierifu,plate,ifu)

    sanitize(gvrot)
    sanitize(gvrote)
    sanitize(svrot)
    sanitize(svrote)
    sanitize(ad2)
    sanitize(ad2e)
    adr = False #set to normalize by gas rotation speed
    ad2 = gvrot**2 - svrot**2

    #make useful arrays
    Re = np.around(np.arange(.1,1.6,.1),1).astype(str)
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified','','','SF Control')
    c = ['k', 'b', 'r', 'gold', 'g', 'm', 'b','','b']

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 4
    control, weights = makecontrol(bpt, lier, c=conc, counts=True,
            mass = lier[masstype])
    bpt[control] = 8

    plt.figure(figsize = (12,4))
    plt.subplot(133)
    
    mass8 = lier[masstype][bpt==8]
    mass2 = lier[masstype][bpt==2]
    if conc == 'pet':
        conc8 = lier[bpt==8]['PETRO_R90']/lier[bpt==8]['PETRO_RE']
        conc2 = lier[bpt==2]['PETRO_R90']/lier[bpt==2]['PETRO_RE']
    elif conc == 'ser':
        conc8 = lier[bpt==8]['SER_N']
        conc2 = lier[bpt==2]['SER_N']

    plt.plot(mass8, conc8, 'b.',
            label = 'control')
    plt.plot(mass2, conc2, 'r.',
            label = 'clier')
    plt.legend()
    plt.title('Control Sample with cLIERS')
    plt.xlabel('Mass')
    plt.ylabel('Concentration')
    #for i in range(len(control)):
    #    plt.text(lier[masstype][control][i], lier['SER_N'][control][i],
    #            str(i))
    #for i in range((bpt==2).sum()):
    #    plt.text(mass2[i], conc2[i],
    #            str(i))

    plt.subplot(131)
    plt.hist(mass8, label  = 'sfcw', histtype = 'step',
            range = (9,12), weights = weights, color = 'b')
    plt.hist(mass8, label  = 'sfcu', histtype = 'step',
            range = (9,12), color = 'cornflowerblue')
    plt.hist(mass2, label = 'cl', histtype = 'step',
            range = (9,12), color = 'r')
    plt.axvline(x = mass8.mean(), color = 'cornflowerblue')
    plt.axvline(x=np.average(mass8, weights = weights), color
            = 'blue')
    plt.axvline(x=mass2.mean(), color = 'r', ls = '--')
    plt.legend()
    plt.title('Mass Histogram')
    plt.xlabel('Mass')
    print('sfw %g sfu %g cl %g' % (np.average(mass8, weights
        = weights), mass8.mean(),
        mass2.mean()))

    plt.subplot(132)
    plt.hist(conc8, label  = 'sfcw', histtype = 'step',
            range = (0,6), weights = weights, color = 'b')
    plt.hist(conc8, label  = 'sfcu', histtype = 'step',
            range = (0,6), color = 'cornflowerblue')
    plt.hist(conc2, label = 'cl', histtype = 'step',
            range = (0,6), color = 'r', ls = '--')
    plt.axvline(x = conc8.mean(), color = 'cornflowerblue')
    plt.axvline(x=np.average(conc8, weights = weights), color
            = 'b')
    plt.axvline(x=conc2.mean(), color = 'r', ls = '--')
    plt.legend()
    plt.xlabel('Concentration')
    plt.title('Concentration Histogram')
    print('sfw %g sfu %g cl %g' % (np.average(conc8, weights
        = weights), conc8.mean(),
        conc2.mean()))
    #print(lier[masstype][bpt==1][:10], mass2[:10])

    plt.tight_layout()
    plt.draw()

    #put the weights in the right shape to deal with
    wa = np.zeros_like(bpt)
    wa[control] = weights
    wahyb = wa[hybtolier]

    #parameters of interest
    interested = [8,2]
    lcut = 11.5
    ucut = 8.5
    median = False
    nbins = 3
    #bins = np.linspace(ucut,lcut,nbins+1)
    bins = makeqbins(nbins, Mi, bpt, hybtolier, interested)
    print(bins)

    #cut the cLIERs into two populations as seen in lierhist
    if 6 in interested:
        r = 3 #radius in Re list to make cut on 
        dip = -.6 #dip in bimodal histogram to cut on
        group = 2 #type number of galaxies to cut (1=sf, 2 = clier)
        ccut(bpt, dip, group, lier, spx, r)

    #more useful arrays
    adprof  = np.zeros((len(Re), len(interested), nbins))
    adeprof = np.zeros((len(Re), len(interested), nbins))
    harcprof  = np.zeros((len(Re), len(interested), nbins))
    harceprof = np.zeros((len(Re), len(interested), nbins))
    strcprof  = np.zeros((len(Re), len(interested), nbins))
    strceprof = np.zeros((len(Re), len(interested), nbins))

    #pick out corresponding Re values
    for j in range(len(Re)):
        ad  = ad2[:,j]
        ade = ad2e[:,j]
        harc  = gvrot[:,j]
        harce = gvrote[:,j]
        strc  = svrot[:,j]
        strce = svrote[:,j]

        if adr:
            ad = ad/(harc**2)
            ade = np.sqrt((ade/(harc**2))**2 + ((2*ad*harce)/(harc**3))**2)
        
        #pick out appropriate bpt type
        for k in range(len(interested)):
            c1 = (bpt==interested[k])[hybtolier]
            print('{0:4} | {1:12}'.format(Re[j], types[interested[k]]),
                    end = '')
            Micut  = Mi[c1]
            Miecut = Mie[c1]
            adcut  = ad[c1]
            adecut = ade[c1]
            harccut  = harc[c1]
            harcecut = harce[c1]
            strccut  = strc[c1]
            strcecut = strce[c1]
            weightscut = wahyb[c1]

            #put data into mass bins
            for l in range(nbins):
                cut = (Micut < bins[l+1]) * (Micut > bins[l])
                print('| {0:3} '.format(cut.sum()), end = '')
                adcutl  = adcut[cut]
                adecutl = adecut[cut]
                Micutl  = Micut[cut]
                Miecutl = Miecut[cut]
                harccutl  = harccut[cut]
                harcecutl = harcecut[cut]
                strccutl  = strccut[cut]
                strcecutl = strcecut[cut]
                weightscutl = weightscut[cut]
        
                if interested[k] == 8:
                    #calculate averages and erros
                    adprof[j,k,l], adeprof[j,k,l] = weighted_avg_std(adcutl, 
                            weightscutl)
                    harcprof[j,k,l], harceprof[j,k,l] = weighted_avg_std(harccutl, 
                            weightscutl)
                    strcprof[j,k,l], strceprof[j,k,l] = weighted_avg_std(strccutl, 
                            weightscutl)
                else:
                    #calculate averages and erros
                    adprof[j,k,l], adeprof[j,k,l] = weighted_avg_std(adcutl, 
                            adecutl)
                    harcprof[j,k,l], harceprof[j,k,l] = weighted_avg_std(harccutl, 
                            harcecutl)
                    strcprof[j,k,l], strceprof[j,k,l] = weighted_avg_std(strccutl, 
                            strcecutl)

                #do median instead
                label = ''
                if median:
                    adprof[j,k,l] = np.nanmedian(adcutl)
                    harcprof[j,k,l] = np.nanmedian(harccutl)
                    strcprof[j,k,l] = np.nanmedian(strccutl)
                    label = 'Median '
            print()
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
    #dist = drp['nsa_zdist'][hybtodrp] * const.c.to('km/s').value / .07
    #hrot = hyb['shrot'] * dist / (206265 * re)
    #hrote = hyb['shrote'] * dist / (206265 * re)
    hrot = hyb['shrot']/re
    hrote = hyb['shrote']/re
    #print(hrot.shape)
    plt.figure()
    for k in range(len(interested)):
        ulim = 2
        rcut = ((bpt==interested[k])[hybtolier] * (hrot < ulim))
        hr = hrot[rcut]
        hre = hrote[rcut]

        plt.hist(hr, bins = 5*ulim, range = (0,ulim), histtype = 'step', 
                label = types[interested[k]], color = c[interested[k]])
        havg, hstd = weighted_avg_std(hr, hre)

        plt.axvline(havg, color = c[interested[k]], ls = '-')
        plt.axvline(havg - hstd, color = c[interested[k]], ls = '--')
        plt.axvline(havg + hstd, color = c[interested[k]], ls = '--')

    plt.title('Stellar Rotation Scale')
    plt.xlabel(r'Rotation Scale (R$_e$)')
    plt.legend()

    #make stuff for plot
    plt.figure(figsize=(12,4))
    c2 = make_cmap(nbins, 'gnuplot2')
    handles = []
    ls = []

    #make lines and labels for legend
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
    if 8 in interested:
        solid = mlines.Line2D([],[], color = 'k', label = 'SF Control')
        handles += [solid]
        ls += ['-']
    
    #go through all types and radii and plot profiles for gas, str, ad
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
                #ax.set_ylim((0, 400))
                ax.tick_params(left = False, labelleft = False)
                plt.grid(True)

    #print(adprof)
    plt.legend(handles = handles)
    plt.tight_layout()
    plt.show()
