#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *

if __name__ == '__main__':
    spx = readall()
    q = (0,25,50,75,100)
    Re = [.25,.50,.75,1.00,1.25]
    cmap = make_cmap(4,'gnuplot2')
    plt.figure(figsize = (12,8))
    interested = [0,1,2,3]

    for i in range(len(Re)):
        plt.subplot(231+i)
        ad = np.log10(spx[i]['ad2_em'])
        ade = spx[i]['ad2_se']
        nmr = spx[i]['elpetro_absmag'][:,1] - spx[i]['elpetro_absmag'][:,4]
        nmre = np.sqrt(spx[i]['elpetro_abmerr'][:,1]**2 + 
                spx[i]['elpetro_abmerr'][:,4]**2)
        mass = np.log10(spx[i]['elpetro_mass'])

        quarts = np.percentile(ad, q)
        for j in interested:
            cut = (ad > quarts[j]) * (ad < quarts[j+1])
            #plt.plot(mass[cut], nmr[cut], '.', c = cmap[j])
            plt.errorbar(mass[cut], nmr[cut], yerr = nmre[cut], fmt = '.', 
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
    plt.show()
