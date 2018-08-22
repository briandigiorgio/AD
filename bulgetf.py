#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from adtools import *

if __name__ == '__main__':
    spx = readall()
    lg = loadfits('prob_table_nb4.fits')
    drp = loadfits('drpall-v2_0_1.fits')

    r = [3]
    for j in r:
        plate = spx[j]['plate']
        ifu = spx[j]['ifudesign']
        drpplate = drp['plate']
        drpifu = [int(i) for i in drp['ifudsgn']]
        
        drpra = drp['objra']
        drpdec = drp['objdec']
        lgra = lg['SDSS_RA']
        lgdec = lg['SDSS_DEC']

        spxtodrp = plateifu(drpplate, drpifu, plate, ifu)
        drptolg,d = match_cats(drpra, drpdec, lgra, lgdec)
        good = (d < .1 * u.arcminute)
        goodi = np.where(good)
        drptolgd = np.delete(drptolg, goodi)

        '''
        plt.plot(lg['Z'][drptolg[good]], drp['nsa_z'][good], '.', alpha = .5)
        plt.xlim((0,.05))
        plt.ylim((0,.15))
        plt.show()
        '''

        ad = spx[j]['ad2_em']
        ade =spx[j]['ad2_se']
        inc = spx[j]['sini']
        ince = spx[j]['sini_e']
        harc = spx[j]['harc_em']
        vrot = spx[j]['gas_vrot']/inc
        harce = spx[j]['harc_se']
        vrote = np.sqrt((spx[j]['gas_vrot_e']/inc)**2 +\
                (spx[j]['gas_vrot']*ince/(inc**2))**2)
        Mi = spx[j]['elpetro_absmag'][:,5]
        Mie= spx[j]['elpetro_abmerr'][:,5]
        mass = np.log10(spx[j]['elpetro_mass'])
        masse = np.zeros_like(Mi)
        bulge = lg['BULGE_TO_TOT_I'][(drptolg * good)[spxtodrp]]
        #plt.hist(bulge, bins = 100, range = (0,1))
        #plt.show()

        #plt.scatter(mass, vrot, s = 10, cmap = 'gnuplot2', vmax = .6, c = bulge)

        bulgeg = np.delete(bulge, np.where(bulge==np.median(bulge)))
        massg = np.delete(mass, np.where(bulge==np.median(bulge)))
        vrotg = np.delete(vrot, np.where(bulge==np.median(bulge)))
        vroteg = np.delete(vrote, np.where(bulge==np.median(bulge)))
        med = np.median(bulgeg)

        for i in range(len(lg)):
            if lg['BULGE_TO_TOT_I'][i] == med:
                print(i,lgra[i], lgdec[i])
        lb = (bulgeg <= med)
        hb = (bulgeg > med)
        print(np.sum(lb))
        print(np.sum(hb))
        plt.errorbar(massg[lb], vrotg[lb], yerr = vroteg[lb], elinewidth = .2,
                ecolor = '.25', fmt = '.', label = 'Low Bulge Fraction',
                color = 'b')
        plt.errorbar(massg[hb], vrotg[hb], yerr = vroteg[hb], elinewidth = .2,
                ecolor = '.25', fmt = '.', label = 'High Bulge Fraction', 
                color = 'orange')

        ax = plt.gca()
        ax.set_yscale('log')
        plt.xlabel('Mass')
        plt.ylabel(r'V$_{rot}$')
        #plt.colorbar(label = 'Bulge Fraction')
        plt.tight_layout()
        plt.legend(title = 'Median = %g' % med)
        plt.show()
