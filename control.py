#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial
from adtools import *

if __name__ == '__main__':
    lier = loadfits('gal_list_v2_0_1_bpt_classify3.fits')
    spx = readall()

    #make useful arrays
    Re = ('0.25', '0.50', '0.75', '1.00', '1.25')
    types = ('Lineless', 'Star-Forming', 'cLIER', 'eLIER', 'AGN',
    'Unclassified', 'cLIER below cut', 'SF below cut')

    #clean up the nans in the bpt designations
    bpt = lier['BPT_C']
    bpt[np.isnan(bpt)] = len(types) - 1
    sf = lier[bpt == 1]
    cl = lier[bpt == 2]

    sfcoords = np.stack((sf['MASS'], sf['SER_N'])).T
    clcoords = np.stack((cl['MASS'], cl['SER_N'])).T

    tree = spatial.KDTree(sfcoords)
    control = tree.query(clcoords)[1]
    bpt[control] = 8

    #print(control)
    #plt.plot(sfcoords[control][:10,0], sfcoords[control][:10,1], 'b.')
    #plt.plot(clcoords[:10,0], clcoords[:10,1], 'r.')
    #plt.show()
