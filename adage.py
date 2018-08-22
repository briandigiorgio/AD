#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from adtools import *

if __name__ == '__main__':
    spx = loadfits('SPX-GAU-MILESHC-composite_1.00Re.fits')
    ff1 = loadfits('manga_firefly-v2_1_2.fits', 1)
    ff2 = loadfits('manga_firefly-v2_1_2.fits', 2)
    drpall = loadfits('drpall-v2_0_1.fits')
    d4 = np.load('wd4000.npy')
    dn4 = np.load('wdn4000.npy')
    lier = loadfits('gal_list_v2_0_1_bpt_classify3.fits')

    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified', 'SF Control')
    c = ['k','b','r','gold','g','m','purple']

    plate = spx['plate']
    ifu = spx['ifudesign']
    ffplate = ff1['PLATE'].astype(int)
    ffifu = np.asarray([int(i) for i in ff1['IFUDSGN']])
    drpplate = drpall['plate']
    drpifu = np.asarray([int(i) for i in drpall['ifudsgn']])
    lierplate = lier['PLATE']
    lierifu = lier['IFUDESIGN']

    spxtoff = plateifu(ffplate, ffifu, plate, ifu)
    spxtodrp = plateifu(drpplate, drpifu, plate, ifu)
    spxtolier = plateifu(lierplate, lierifu, plate, ifu)
    
    ad = spx['ad2_em']
    ade = spx['ad2_se']
    #age = ff2['LW_Z_1Re'][spxtoff]
    agee = ff2['LW_AGE_1Re_ERROR'][spxtoff]
    dn4000 = dn4[2]
    sern = drpall[spxtodrp]['nsa_sersic_n']
    harc = lier['EWHA_RE'][spxtolier]
    age = sern

    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = 5
    #control = makecontrol(bpt,lier)
    #bpt[control] = 6

    #mass = spx['elpetro_mass']
    #ffmass = ff1['PHOTOMETRIC_MASS'][spxtoff]
    #plt.semilogx(mass, ffmass, 'b.')

    #plt.semilogy(d4000, ad, 'b.')
    #plt.errorbar(age, ad, xerr = agee, yerr = ade, fmt = '.', ecolor = '.75',
    #        elinewidth = .2)

    plt.figure(figsize = (10,4))
    plt.subplot(121)
    for i in range(len(ad)):
        plt.plot(dn4000[i],ad[i],'.',c=c[bpt[spxtolier][i].astype(int)])
    handles = []
    for k in range(len(types)):
        handles += [mlines.Line2D([], [], color = c[k], label = types[k], 
                marker = '.')]
    handles = handles[:4] + [handles[5]]
    plt.legend(handles = handles)
    plt.xlabel(r'D$_n$4000')
    plt.ylabel(r'AD$^2$')
    ax = plt.gca()
    ax.set_yscale('log')
    plt.grid(True)

    plt.subplot(122)
    plt.scatter(dn4000, ad, c=age, cmap='gnuplot', s=10)#, vmin = -.5)#, vmax = 0)
    plt.colorbar(label = 'Sersic n')
    plt.xlabel(r'D$_n$4000')
    ax = plt.gca()
    ax.set_yscale('log')

    plt.grid(True)
    plt.tight_layout()
    plt.show()
