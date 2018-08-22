#!/usr/bin/env python

import numpy
from astropy.io import fits
from adtools import *

def get_bot(plate, ifu):
    # plate is the plate number for all the galaxies of interest
    # ifu is the ifudesign number for all the galaxies of interest

    # Read the matching catalog
    matchfits = 'manga_catalog_match.fits.gz'
    match_hdu = fits.open(matchfits)
    ndrp = len(match_hdu['DRPMATCH'].data['PLATE'])
    # In the DRPMATCH extension, there is one table row per DRPall entry

    # This gets the indices of the rows in the DRPall file with the
    # selected plate and ifu
    indx=numpy.array([numpy.arange(ndrp)[(match_hdu['DRPMATCH'].data['PLATE']
        ==p) & (match_hdu['DRPMATCH'].data['IFUDESIGN'] == f)][0]
                                                for p,f in zip(plate, ifu)])

    # Of the selected DRPall rows, find the ones with a matching
    # measurement in the Lackner & Gunn catalog; LG_INDX is the index of
    # the row in the 'lackner_gunn_2012.fits.gz' binary tables; any
    # value that =-1 does *not* have a match
    lg_indx = match_hdu['DRPMATCH'].data['LG_INDX'][indx] > -1
    print(indx.shape, lg_indx.shape, match_hdu['DRPMATCH'].data['PLATE'].shape)

    # Read the Lackner & Gunn data
    lgfits = 'lackner_gunn_2012.fits.gz'
    lghdu = fits.open(lgfits)

    # For each input plate-ifu, get the LG-selected model type
    lg_fit_selection = numpy.ma.masked_all(len(plate), dtype='<U3')
    lg_fit_selection[lg_indx] \
        = lghdu['SUMMARY'].data['TYPE'][match_hdu['DRPMATCH'].data['LG_INDX'][indx[lg_indx]]]

    # Get the bulge-disk decomposition using a de Vaucouleurs bulge
    lg_nb4_bot = numpy.ma.masked_all(len(plate), dtype=float)
    lg_nb4_bot[lg_indx] \
        = lghdu['NB4'].data['BULGE_TO_TOT'][match_hdu['DRPMATCH'].data['LG_INDX'][indx[lg_indx]],2]
    lg_nb4_bot[lg_nb4_bot < 0] = numpy.ma.masked
    
    # Get the bulge-disk decomposition using an exponential bulge
    lg_nb1_bot = numpy.ma.masked_all(len(plate), dtype=float)
    lg_nb1_bot[lg_indx] \
        = lghdu['NB1'].data['BULGE_TO_TOT'][match_hdu['DRPMATCH'].data['LG_INDX'][indx[lg_indx]],2]
    lg_nb1_bot[lg_nb1_bot < 0] = numpy.ma.masked

    # Return the selection and the B/T for each fit; the length of the
    # returned arrays is the same as the input plate and ifu lists; any
    # plate-ifu observation not matched to an LG measurement is masked
    return lg_fit_selection, lg_nb4_bot, lg_nb1_bot

if __name__ == "__main__":
    spx = loadfits('SPX-GAU-MILESHC-composite_1.00Re.fits')
    get_bot(spx['plate'], spx['ifudesign'])
