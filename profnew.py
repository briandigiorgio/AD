#!/usr/bin/env python

from adtools import *

def was(values, errors):
    avg = np.average(values, weights = 1/errors**2, axis = 0)
    err = np.sqrt(np.average((values - avg)**2, axis = 0, 
        weights = 1/(errors**2)))
    return avg, err

def waso(values, errors):
    avg = np.average(values, weights = 1/errors, axis = 0)
    var = np.sum(errors * (values - avg)**2)/errors.sum()
    return avg, np.sqrt(var/len(values))

if __name__ == '__main__':
    hyb = loadfits('adsample.fits')
    hyb = hyb[hyb['cut']]
    lier = loadfits('MPL-6_master_catalogue_Apr9.fits')
    bpt = lier['BPT_C']
    hybtolier = np.load('hybtodrplier.npy')[1]
    masstype = 'MASS_SER'#'MASS_ELL_PETRO'
    conc = 'ser'#'pet'

    svrot = np.ma.masked_array(hyb['svrot'], hyb['svrot'] < -100)#== -9999)
    svrote = np.ma.masked_array(hyb['svrote'], hyb['svrote'] < -100)#== -9999)
    gvrot = np.ma.masked_array(hyb['gvrot'], hyb['gvrot'] < -100)#== -9999)
    gvrote = np.ma.masked_array(hyb['gvrote'], hyb['gvrote'] < -100)#== -9999)
    ad = np.sqrt(gvrot**2 - svrot**2)
    ade = np.sqrt(((gvrot*gvrote)**2 + (svrot*svrote)**2)/(gvrot**2 - svrot**2))
    Mi = lier[masstype][hybtolier]

    #control, weights = makecontrol(bpt, lier, c = 'pet', counts = True, stl = hybtolier)
    control, weights = make_control(bpt, lier, masstype, conc)
    wa = np.zeros_like(bpt)
    wa[control] = weights
    wa = wa[hybtolier]
    bpt[control] = 8

    interested = [2,8]
    nbins = 3
    bins = makeqbins(nbins, Mi, bpt, hybtolier, interested)
    c = make_cmap(nbins, 'gnuplot2', invert = True)
    Re = np.arange(.1,1.6,.1).round(1)

    ls = ['--','-']
    labels = ['cLIER', 'SF Control']
    arrays = [gvrot, svrot, ad]
    errors = [gvrote, svrote, ade]
    titles = ['gvrot', 'svrot', 'ad']

    plt.figure(figsize = (4*len(arrays),4))
    for j in range(len(arrays)):
        plt.subplot(int('1%s%s' % (len(arrays), j+1)))
        plt.title(titles[j])
        plt.xlabel(r'R$_e$')
        plt.ylabel(titles[j])
        for i in range(nbins):
            for k in range(len(interested)):
                cut = bpt[hybtolier] == interested[k] * (Mi > bins[i]) \
                        * (Mi <= bins[i+1])
                print(labels[k], titles[j], bins[i].round(2), bins[i+1].round(2), 
                        cut.sum())
                #print(np.average(arrays[j][cut], axis = 0, weights =
                #    1/(errors[j][cut])**2))
                #if j==0 and i==2:
                #    d = arrays[j][cut][:,2]
                #    e = errors[j][cut][:,2]
                #    print(d[d<0], e[d<0])
                #    plt.hist(d)
                #    plt.show()
                #avg, err = waso(arrays[j][cut], errors[j][cut])
                avg = arrays[j][cut].mean(axis = 0)
                
                plt.plot(Re, avg, c = c[i], ls = ls[k], label = labels[k])
                #plt.fill_between(Re, avg-err, avg+err, color = c[i], alpha = .1)
        plt.legend()
    plt.tight_layout()
    plt.show()
