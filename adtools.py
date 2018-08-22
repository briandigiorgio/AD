import time
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors, rc, cm
import argparse
import scipy.stats as stats
from scipy.optimize import curve_fit
from astropy import coordinates
import astropy.units as u
from scipy import spatial
from urllib.request import urlretrieve
from PIL import Image

def load_fits(file, i = 1):
    #return fits.open(file)[i].data.copy()
    hdu = fits.open(file)
    data = hdu[i].data.copy()
    hdu.close()
    return data

def loadfits(file, i=1):
    return load_fits(file, i)

def lex(file, i=1):
    hdu = loadfits(file, i=i)
    print('\n'.join(hdu.columns.names))
    return hdu

#reads fits file of data from Kyle, spits out relevant data arrays
#should be of the form 'SPX-GAU-MILESHC-composite_*.fits'
#largely borrowed from Kyle
def read_ad_mass(file, adv = False):

    #loads file
    data = load_fits(file)

    #names of columns to be pulled
    mag = 'elpetro_absmag'
    mage = 'elpetro_abmerr'

    #pull Mi data
    Mi = data[mag][:,5]
    Mie = data[mage][:,5]
    
    #divide by gas rotation speed if user specified
    if adv:
        ad = data['ad2_em']/data['harc_em']
        ade = np.sqrt((data['ad2_se']/data['harc_em'])**2 +
              ((data['ad2_se']*data['harc_se'])/(data['harc_em']**2))**2)

    else: 
        ad = data['ad2_em']
        ade = data['ad2_se']
    
    return Mi, Mie, ad, ade

#reads info on magnitude, ad, color, rotation speed, and identifying info
#taken mostly from Kyle, works on SPX... files
def read_ad_mag(file):

    #loads file
    data = load_fits(file)

    #names of columns to be pulled
    mag = 'elpetro_absmag'
    mage = 'elpetro_abmerr'
    mass = 'elpetro_mass'
    
    #pull Mi data
    Mi = data[mag][:,5]
    Mie = data[mage][:,5]

    #calculate N-r
    Nmr = data[mag][:,1] - data[mag][:,4]
    Nmre = np.sqrt(np.square(data[mage][:,1]) + np.square(data[mage][:,4]))

    #pull rotational velocity data
    gasrc = data['harc_em']
    gasrce = data['harc_se']

    ad = data['ad2_em']
    ade = data['ad2_se']

    plate = data['plate'].astype(str)
    ifu = data['ifudesign'].astype(str)

    return plate, ifu, Mi, Mie, Nmr, Nmre, gasrc, gasrce, ad, ade

def find_repeats(array):
    return np.setdiff1d(array, np.unique(array))

def powerlaw(x, a, b, c):
    return a * np.power(x, b) + c

def line(x, m, b):
    return m * x + b

def exponential(x, a, b):
    return a * b ** x

def weighted_avg_std(values, errors):
    if len(values) == 0 or len(errors) == 0 or np.sum(errors) == 0:
        return (np.nan, np.nan) #failure

    if np.isnan(values).any():
        values = np.ma.array(values, mask = np.isnan(values))
        errors = np.ma.array(errors, mask = np.isnan(errors))

    average = np.average(values, weights=1/errors)
    sumerrors = np.sum(errors)
    variance = np.sum((errors * (values - average)**2))/sumerrors
    return (average, np.sqrt(variance)/np.sqrt(len(values)))

#makes an array of length numlines going through the specified color map
#plt.plot(..., c = colors[i]) when plotting multiple lines
def make_cmap(numlines, cmap, invert = False):
    cnorm = colors.Normalize(vmin = 0, vmax = numlines)
    scalarmap = cm.ScalarMappable(norm = cnorm, cmap = cmap)
    if invert:
        return scalarmap.to_rgba(range(numlines)[::-1])
    else:
        return scalarmap.to_rgba(range(numlines))

#returns an array of indices of cat2 that correspond to the same item in cat1
def match_cats(ra1, dec1, ra2, dec2):
    cat1 = coordinates.SkyCoord(ra = ra1 * u.degree, dec = dec1 * u.degree)
    cat2 = coordinates.SkyCoord(ra = ra2 * u.degree, dec = dec2 * u.degree)
    #cat1tocat2, d2d, d3d = cat1.match_to_catalog_sky(cat2)
    cat1tocat2, d2d, d3d = coordinates.match_coordinates_sky(cat1, cat2)
    '''
    print(len(cat1tocat2))
    bad = []
    for i in range(len(d2d)):
        if d2d[i].is_within_bounds('10m',None):
            bad += [i]
    cat1tocat2 = np.delete(cat1tocat2, bad)
    print(len(cat1tocat2), len(bad))
    '''
    return cat1tocat2, d2d

#reads fits file of data from Kyle, spits out relevant data arrays
#should be of the form 'SPX-GAU-MILESHC-composite_*.fits'
#largely borrowed from Kyle
def read_mass(file, adv = False):

    #loads file
    data = load_fits(file)

    #names of columns to be pulled
    mass = 'elpetro_mass'
    #mass = 'sersic_mass'

    #pull Mi data
    m = data[mass]
    
    #divide by gas rotation speed if user specified
    if adv:
        ad = data['ad2_em']/data['harc_em']
        ade = np.sqrt((data['ad2_se']/data['harc_em'])**2 +
              ((data['ad2_se']*data['harc_se'])/(data['harc_em']**2))**2)

    else: 
        ad = data['ad2_em']
        ade = data['ad2_se']
    
    return m, ad, ade

#loads all of the spx files for different radii and puts them in an array
def readall(f25 = 'SPX-GAU-MILESHC-composite_0.25Re.fits',
        f50 = 'SPX-GAU-MILESHC-composite_0.50Re.fits', 
        f75 = 'SPX-GAU-MILESHC-composite_0.75Re.fits', 
        f10 = 'SPX-GAU-MILESHC-composite_1.00Re.fits', 
        f12 = 'SPX-GAU-MILESHC-composite_1.25Re.fits', others = []):
    spx25 = loadfits(f25)
    spx50 = loadfits(f50)
    spx75 = loadfits(f75)
    spx10 = loadfits(f10)
    spx12 = loadfits(f12)
    return np.asarray([spx25, spx50, spx75, spx10, spx12]+others)

#make a cut in the ad/hrot ratio between the two populations from lierhist
#bpt = array from francesco, cut = where to split the populations
#group = group in bpt to split, spx = array of all Re files, 
#r = which Re to use for the cut, l and u for a certain mag range
def ccut(bpt, dip, group, lier, spx, r, i=6, l=0, u=0):
    plate = spx[r]['plate'].astype(str)
    ifu = spx[r]['ifudesign'].astype(str)
    plateifuspx = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])

    lierplate = lier['PLATE'].astype(str)
    lierifu = lier['IFUDESIGN'].astype(str)
    plateifulier = np.asarray([lierplate[i] + lierifu[i] 
        for i in range(len(lierplate))])

    spxtolier = np.zeros(len(plateifuspx))
    for l in range(len(plateifuspx)):
        for m in range(len(plateifulier)):
            if plateifuspx[l] == plateifulier[m]:
                spxtolier[l] = m
    spxtolier = spxtolier.astype(int)

    ad10 = spx[r]['ad2_em']
    harc10  = spx[r]['harc_em']
    cut = spxtolier * (bpt[spxtolier] == group) * \
            (np.log10(ad10/(harc10**2)) < dip)
    
    if l:
        print(np.sum(spx[r]['elpetro_absmag'][:,5] > l))
        cut = cut * (spx[r]['elpetro_absmag'][:,5] < l)
    if u:
        print(np.sum(spx[r]['elpetro_absmag'][:,5] < u))
        cut = cut * (spx[r]['elpetro_absmag'][:,5] > u)

    bpt[cut] = i

def ccutul(bpt, dip, group, lier, spx, r, l, u, i=6):
    plate = spx[r]['plate'].astype(str)
    ifu = spx[r]['ifudesign'].astype(str)
    plateifuspx = np.asarray([plate[i]+ifu[i] for i in range(len(plate))])

    lierplate = lier['PLATE'].astype(str)
    lierifu = lier['IFUDESIGN'].astype(str)
    plateifulier = np.asarray([lierplate[i] + lierifu[i] 
        for i in range(len(lierplate))])

    spxtolier = np.zeros(len(plateifuspx))
    for l in range(len(plateifuspx)):
        for m in range(len(plateifulier)):
            if plateifuspx[l] == plateifulier[m]:
                spxtolier[l] = m
    spxtolier = spxtolier.astype(int)
    cutl = (spx[r]['elpetro_absmag'][:,5] < -21)
    cutu = (spx[r]['elpetro_absmag'][:,5] > -23)

    ad10 = spx[r]['ad2_em']
    harc10  = spx[r]['harc_em']
    cut = spxtolier * (bpt[spxtolier] == group) * \
            (np.log10(ad10/(harc10**2)) < dip) * cutl * cutu
    print(np.sum(cut.astype(bool)))
    bpt[cut] = i

#makes a colorbar that is also a histogram
def colorhist(data, cmap, ax = None, bins = 10, range = None, density = False, label = None):
    if ax:
        n,b,p = ax.hist(data, bins = bins, range = range, density = density,
                label = label)
    else:
        n,b,p = plt.hist(data, bins = bins, range = range, density =
                density, label = label)

    bincenters = .5 * (b[:-1] + b[1:])
    col = bincenters - np.min(bincenters)
    col /= np.max(col)
    cmap = plt.cm.get_cmap(cmap)


    for ci, pi in zip(col, p):
        plt.setp(pi, 'facecolor', cmap(ci))

    return n,b,p    

#makes a control sample of SF galaxies that are similar to cLIERs
#uses a mass parameter and a concentration parameter
def makecontrol(bpt, lier, c='pet', mass = np.asarray([0]), stl =
        np.asarray([0]), spx = np.asarray([0]), norm = False, counts = False,
        allbpt = False):

    #if fed spxtolier array, reorder things so mass can be used
    if stl.any():
        lier = lier[stl]
        bpt = bpt[stl]
        sf = lier[bpt==1]
        cl = lier[bpt==2]
        if mass.any():
            mass = mass[stl]
    else:
        #load starforming and cLIER galaxies
        sf = lier[bpt == 1]
        cl = lier[bpt == 2]

    #make a range of indices for cutting
    indx = np.arange(len(lier))

    #get the appropriate concentration parameter (sersic or petrosian)
    if c == 'ser':
        concsf = sf['SER_N']
        conccl = cl['SER_N']
    if c == 'pet':
        #try two different things for different MPL5 and 6 versions
        try:
            concsf = sf['PETRO_90']/sf['PETRO_50']
            conccl = cl['PETRO_90']/cl['PETRO_50']
        except:
            concsf = sf['PETRO_R90']/sf['PETRO_RE']
            conccl = cl['PETRO_R90']/cl['PETRO_RE']

    #get desired mass numbers for the galaxies
    if mass.any():
        masssf = mass[bpt==1]
        masscl = mass[bpt==2]
    else:
        try:
            masssf = sf['MASS']
            masscl = cl['MASS']
        except:
            if c == 'pet':
                masssf = sf['MASS_ELL_PETRO']
                masscl = cl['MASS_ELL_PETRO']
            elif c == 'ser':
                masssf = sf['MASS_SER']
                masscl = cl['MASS_SER']

    sfgood = np.isfinite(masssf)# * np.isfinite(concsf)
    clgood = np.isfinite(masscl)# * np.isfinite(conccl)

    #normalize the coordinates so the matching isn't blown out of the water by 
    #gigantic mass numbers
    if norm:
        #sfgood = np.ones_like(masssf, dtype = bool)
        #clgood = np.ones_like(masscl, dtype = bool)
        #print(np.sum(np.isnan(masssf)))
        masssf = normalize(masssf[sfgood])
        masscl = normalize(masscl[clgood])
        concsf = normalize(concsf[sfgood])
        conccl = normalize(conccl[clgood])
    else:
        masssf = masssf[sfgood]
        masscl = masscl[clgood]
        concsf = concsf[sfgood]
        conccl = conccl[clgood]


    #make ordered pairs
    sfcoords = np.stack((masssf, concsf)).T
    clcoords = np.stack((masscl, conccl)).T

    #find the closest SF to each cLIER to get a control
    control = spatial.KDTree(sfcoords).query(clcoords)[1]

    #plt.figure(figsize = (12,4))
    #plt.subplot(133)
    #plt.plot(masssf[control][:], lier[bpt == 1][control]['SER_N'][:], 
    #        'b.')
    #plt.plot(masscl[:], lier[bpt == 2]['SER_N'][:],'r.')
    #for i in range(len(control)):
    #    plt.text(masssf[control][i], lier['SER_N'][bpt==1][control][i],
    #            str(i))
    #for i in range(len(masscl)):
    #    plt.text(masscl[i], lier['SER_N'][bpt==2][i],
    #            str(i))


    if counts:
        crange = np.arange(np.amax(control+1))
        freq = np.bincount(control)
        countdict = dict(zip(crange, freq))
        weights = np.array([countdict[i] for i in control])
        unq = np.unique(control, return_index = True)[1]

        #plt.subplot(131)
        #plt.hist(masssf[control], label  = 'sfcw', histtype = 'step',
        #        weights = weights)
        #plt.hist(masssf[control], label  = 'sfcu', histtype = 'step')
        #plt.hist(masscl, label = 'cl', histtype = 'step')
        #plt.axvline(x = masssf[control].mean(), color = 'orange')
        #plt.axvline(x=np.average(masssf[control], weights = weights), color
        #        = 'b')
        #plt.axvline(x=masscl.mean(), color = 'g')
        #print('sfu %s sfw %s cl %s' % (masssf[control].mean(), np.average(masssf[control], weights =
        #    weights), masscl.mean()))
        #plt.legend()

        #plt.subplot(132)
        #plt.hist(concsf[control], label  = 'sfcw', histtype = 'step',
        #        weights = weights)
        #plt.hist(concsf[control], label  = 'sfcu', histtype = 'step')
        #plt.hist(conccl, label = 'cl', histtype = 'step')
        #plt.axvline(x = concsf[control].mean(), color = 'orange')
        #plt.axvline(x=np.average(concsf[control], weights = weights), color
        #        = 'blue')
        #plt.axvline(x=conccl.mean(), color = 'g', ls = '--')
        #print('sfu %s sfw %s cl %s' % (concsf[control].mean(), np.average(concsf[control], weights =
        #    weights), conccl.mean()))
        #plt.legend()

        ##print(masssf[control][:10], masscl[:10])
        #plt.tight_layout()
        #plt.show()
        return indx[bpt == 1][control][unq], weights[unq]

    return indx[bpt == 1][control]

def make_control(bpt, lier, masstype, conc):
    sf = lier[bpt==1][np.isfinite(lier[bpt==1][masstype])]
    cl = lier[bpt==2][np.isfinite(lier[bpt==2][masstype])]
    masssf = sf[masstype]
    masscl = cl[masstype]
    
    if conc=='ser':
        concsf = sf['SER_N']
        conccl = cl['SER_N']
    elif conc=='pet':
        concsf = sf['PETRO_R90']/sf['PETRO_RE']
        conccl = cl['PETRO_R90']/cl['PETRO_RE']

    sfcoords = np.stack((masssf, concsf)).T
    clcoords = np.stack((masscl, conccl)).T
    control = spatial.KDTree(sfcoords).query(clcoords)[1]

    counts = dict(zip(np.arange((control+1).max()), np.bincount(control)))
    weights = np.array([counts[i] for i in control])
    unique = np.unique(control, return_index = True)[1]

    return control[unique], weights[unique]

    




#matches the plate/ifu indices of catalogs (1=bigger, 2=smaller)
#I guess this should work for any two int/float coordinate types
def plateifu(plate1, ifu1, plate2, ifu2):
    plateifubig = np.stack((plate1,ifu1)).T
    plateifusmall = np.stack((plate2,ifu2)).T
    return spatial.KDTree(plateifubig).query(plateifusmall)[1]

def normalize(array):
    return array/np.nanmax(array)

#pull a thumbnail image for an manga galaxy from sdss
#square with side npix, scale is zoom
#makes annoying out file that is unavoidable as far as I can tell
def getimage(ra, dec, scale, npix, out = 'out'):
    url = ("http://skyservice.pha.jhu.edu/DR8/ImgCutout/"
           "getjpeg.aspx?ra=%.8f&dec=%.8f&scale=%.2f&width=%i&height=%i"
           % (ra, dec, scale, npix, npix))
    urlretrieve(url, out)
    return Image.open(out)

#get the ra and dec of a given plate ifu, drpall in loadfits format
def getradec(plate, ifu, drpall):
    drpifu = [int(i) for i in drpall['ifudsgn']]
    todrp = plateifu(drpall['plate'], drpifu, plate, ifu)
    return drpall['objra'][todrp], drpall['objdec'][todrp]

#makes a collage of thumbnails from sdss of a given set of plateifus
#largely lifted from Alexie's plot_utils.py from astro 214
def sdsscollage(plate, ifu, nrow, ncol, npix, scale, drpall):
    ra, dec = getradec(plate, ifu, drpall)
    fig, axs = plt.subplots(nrow, ncol, figsize = (2*ncol, 2*nrow))
    for i in range(len(plate)):
        ax = axs.flatten()[i]
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.imshow(getimage(ra[i], dec[i], scale, npix), origin = 'lower')
        ax.text(.5*npix, .1*npix, '%s-%s' % (plate[i], ifu[i]), color = 'w', 
                fontsize = 12, horizontalalignment = 'center')
    plt.subplots_adjust(hspace = 0, wspace = 0)
    plt.show()

#replaces bad values with good values in an array
#I'm pretty sure there's alread a function for this in numpy
def sanitize(array, value = -9999, replace = np.nan):
    array[array==value] = replace
    return array

#take the svrot/gvrot arrays of uneven arrays and truncates them to just be
#the last 15 items that actually matter and put them into a workable 2d array
def fixnpy(npy, start = -15):
    array = np.load(npy)
    for i in range(len(array)):
        try:
            len(array[i])
        except:
            array[i] = np.zeros_like(array[0])
        array[i] = array[i][start:]
    return np.stack(array)

#make bins by mass quantile, returns boundaries
def makeqbins(nbins, Mi, bpt, hybtolier, interested):
    bincut = np.zeros_like(Mi, dtype = bool)
    for k in range(len(interested)):
        bincut += (bpt == interested[k])[hybtolier].astype(bool)
    #bincut = np.array([(bpt==i)[hybtolier] for i in interested], dtype = bool)
    return np.percentile(Mi[bincut], np.linspace(0,100,nbins+1))
