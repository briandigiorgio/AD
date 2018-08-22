#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from marvin.tools.maps import Maps
from marvin.utils.general.general import downloadList
from marvin import config
from adtools import *
import time


if __name__ == '__main__':
    start = time.time()
    config.setRelease('MPL-5')
    spx = readall()
    #spx = loadfits('SPX-GAU-MILESHC-coposite_1.00Re.fits')
    nbins = 5
    stepsize = .25
    bins = np.arange(0, stepsize*(nbins+1), stepsize)
    alldata = []

    #for i in range(len(spx)):
    for i in [3]:
        print(len(spx[i]))
        plateifu = []
        d4000 = np.zeros(len(spx[i]))
        res = np.zeros_like(d4000)

        for j in range(len(spx[i])):
            plateifu+=[str(spx[i]['plate'][j])+'-'+str(spx[i]['ifudesign'][j])]
        #downloadList(plateifu, dltype = 'map')
            maps = Maps(plateifu = plateifu[j])
            alld4 = maps['specindex-dn4000']
            radius = maps['spx_ellcoo_elliptical_radius'].value
            re = float(maps.header['REFF'])
            res[j] = re
            radius /= re
            
            #for k in range(nbins):
            #    d4000[j,k] = np.average(alld4[(radius < bins[k+1]) * \
            #        (radius > bins[k])])
            #d4000[j]=np.average(alld4[(radius < bins[i+1])*(radius > bins[i])])
            d4000[j]=np.average(alld4[(radius < 1.25)*(radius > .75)])
            
            #alldata += [d4000[j]]
            print(j, plateifu[j], re, d4000[j])

        print(spx[i]['plate'].shape, spx[i]['ifudesign'].shape, d4000.shape)
        alldata = np.stack((spx[i]['plate'], spx[i]['ifudesign'], d4000))
        alldata = np.asarray(alldata)
        plateifu = np.asarray(plateifu)
        np.save('wdn4000',alldata)
        np.save('plateifu',plateifu)
        redata = np.stack((spx[i]['plate'], spx[i]['ifudesign'], res))
        np.save('re',redata)
        print(time.time()-start)
