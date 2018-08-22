#!/usr/bin/env python

import matplotlib as plt
import numpy as np
from adtools import *

if __name__ == '__main__':
    spx = readall()
    Re = [3]

    for j in Re: 
        plate = spx[j]['plate']
        ifu = spx[j]['ifudesign']
        ad = spx[j]['ad2_em']
        ade =spx[j]['ad2_se']
        inc = spx[j]['sini']
        ince = spx[j]['sini_e']
        harc = spx[j]['harc_em']
        vrot = spx[j]['gas_vrot']/inc
        harce = spx[j]['harc_se']
        vrote = np.sqrt((spx[j]['gas_vrot_e']/inc)**2 +\
                (spx[j]['gas_vrot']*ince/(inc**2))**2)
        mass = np.log10(spx[j]['elpetro_mass'])
        vrotm = vrot -  harc
        vrotme = np.sqrt((spx[j]['gas_vrot_e']/inc)**2 +\
                (spx[j]['gas_vrot']*ince/(inc**2))**2 + harce**2)
        re = np.load('re.npy')[2]


        '''
        plt.errorbar(vrot, harc, xerr = vrote, yerr = harce, fmt = '.', 
                elinewidth = .2, ecolor = '.25', ms = 5)
        x = np.linspace(0,500,500)
        plt.plot(x,x,'r--',zorder = 10)

        ax = plt.gca()
        ax.set_xlim((0,450))
        ax.set_ylim((0,450))
        plt.xlabel('gas_vrot/sini')
        plt.ylabel('harc_em')
        ax.set_aspect(1)
        plt.tight_layout()
        plt.show()

        plt.errorbar(mass, vrotm, yerr = vrotme, fmt = '.', elinewidth = .2,
                ecolor = '.25', ms = 5)
        plt.axhline(y=0, color= 'k')
        plt.xlabel('Mass')
        plt.ylabel('gas_vrot/sini - harc_em')
        ax = plt.gca()
        #ax.set_yscale('log', nonposy='clip')
        ax.set_ylim((-50,100))
        plt.tight_layout()
        plt.show()
        '''

        plt.plot(spx[3]['ifu']+2, spx[3]['gas_hrot'], 'b.', ms=1, label =
                'gas_hrot')
        plt.plot(spx[3]['ifu']-2, re, 'g.', ms=1, label = r'R$_e$')
        ifun = [19,37,61,91,127]
        ifusize = np.asarray([12,17,22,27,32])/2
        for i in ifun:
            plt.plot(i+2, np.median(spx[3]['gas_hrot'][spx[3]['ifu'] == i]),
                'r_', ms = 20)
            plt.plot(i-2, np.median(re[spx[3]['ifu'] == i]),
                'm_', ms = 20)
        plt.plot(ifun,ifusize, 'k_', ms = 30)
        plt.xticks(ifun,ifun)
        plt.xlabel('IFU size')
        plt.ylabel('Arcsec')
        plt.legend()
        plt.show()
