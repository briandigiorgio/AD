#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from adtools import *
from astropy.coordinates import Angle
import time

if __name__ == '__main__':
    start = time.time()
    drpall = load_fits('drpall-v2_3_1.fits')
    gz = load_fits('GalaxyZoo1_DR_table2.fits')
    plate,ifu,Mi,Mie,Nmr,Nmre,gasrc,gasrce,ad,ade = \
        read_ad_mag('SPX-GAU-MILESHC-composite_0.25Re.fits')

    #gzra =  Angle(gz['RA'] + ' hours').degree
    #gzdec = Angle(gz['DEC'] + ' deg').degree
    #np.save('gzradec', np.asarray((gzra, gzdec)))
    gzra, gzdec = np.load('gzradec.npy')
    drptogz, d2d = match_cats(drpall['ifura'], drpall['ifudec'], gzra, gzdec)
    gztodrp, d2d2 = match_cats(gzra, gzdec, drpall['ifura'], drpall['ifudec'])
    good = (d2d.arcsecond < 5)
    drpall = drpall[good]
    drptogz = drptogz[good]

    plateifuspx = np.asarray([plate[i] + ifu[i] for i in range(len(plate))])
    drpplate = drpall['plate'].astype(str)
    drpifu = drpall['ifudsgn'].astype(str)
    plateifudrp = np.asarray([drpplate[i]+drpifu[i] 
                  for i in range(len(drpplate))])

    spxtodrp = np.asarray([np.argmax(np.abs(plateifudrp==plateifuspx[i]))
                    for i in range(len(plateifuspx))])
    drptospx = np.asarray([np.argmax(np.abs(plateifudrp[i]==plateifuspx))
                    for i in range(len(plateifudrp))])

    print(len(gz), len(drpall), len(ad))
    print(len(drptogz), len(gztodrp), len(spxtodrp), len(drptospx))
    print(spxtodrp[:20])
    pcs = gz['P_CS_DEBIASED'][drptogz[spxtodrp]]
    pcs2 = gz['P_CS_DEBIASED'][drptogz][spxtodrp]
    #pcs3 = ad[spxtodrp][drptogz]
    plt.semilogy(pcs, ad, 'b,')
    plt.semilogy(pcs2, ad, 'r,')
    #plt.semilogy(g
    print("Time taken: ", time.time()-start)
    plt.show()
