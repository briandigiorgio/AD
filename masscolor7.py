#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *

if __name__ == '__main__':
    cut = np.load('modelfitting/bigcut.npy')
    hyb = np.load('modelfitting/bth__.npy')[cut]
    drpall = loadfits('drpall-v2_4_3.fits')
    ngal = len(cut)
    
    plate = np.load('modelfitting/platescut.npy')
    ifu = np.load('modelfitting/ifuscut.npy')
    drpplate = drpall['plate']
    drpifu = [int(i) for i in drpall['ifudsgn']]
    hybtodrp = plateifu(drpplate,drpifu,plate,ifu)

    ad2 = fixnpy('modelfitting/ad2.npy')
    mass = np.log10(drpall['nsa_elpetro_mass'][hybtodrp])
    nmr = drpall['nsa_elpetro_absmag'][:,1][hybtodrp] -\
            drpall['nsa_elpetro_absmag'][:,4][hybtodrp]
    #nmre = np.sqrt(drpall['nsa_ellpetro_abmerr'][:,1][hybtolier]**2 +\
    #            drpall['nsa_ellpetro_abmerr'][:,4][hybtolier]**2)
    nmre = np.zeros_like(nmr)

    q = (0,25,50,75,100)
    Re = np.around(np.arange(.1,1.6,.3),1)
    cmap = make_cmap(4,'gnuplot2')
    interested = [0,1,2,3]
    interestedi = np.array(interested)[:,np.newaxis]

    #for interested in interestedi:
    plt.figure(figsize = (12,8))
    for i in range(len(Re)):
        plt.subplot(2,3,i+1)
        good = (ad2[:,i] > 0)
        ad = ad2[:,i][good]
        quarts = np.percentile(ad, q)
        print(np.sum(good), quarts)
        #plt.hist(ad)
        #[plt.axvline(x=qi) for qi in quarts]
        #plt.show()
        for j in interested:
            cut = (ad > quarts[j]) * (ad < quarts[j+1])
            #plt.plot(mass[cut], nmr[cut], '.', c = cmap[j])
            plt.errorbar(mass[good][cut], nmr[good][cut], yerr =
                    nmre[good][cut], fmt = '.', 
                    elinewidth = .2, color = cmap[j], ecolor = '.05', 
                    label = '%g - %g%%' % (q[j],q[j+1]))

        #plt.scatter(mass, nmr, c = ad, cmap = 'gnuplot2', s=10)
        #plt.colorbar()
        ax = plt.gca()
        ax.set_xlim((8.5,11.5))
        ax.set_ylim((1,7))
        plt.xlabel('Mass')
        plt.ylabel(r'N-$r$')
        plt.title(r'%g R$_e$' % Re[i])
        plt.legend(title='AD Quartile')
        plt.tight_layout()
        plt.subplots_adjust(hspace=0, wspace = 0)
    plt.show()
