# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

# http://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html
# http://docs.astropy.org/en/stable/io/fits/

"""
=======================================
Read and plot an image from a FITS file
=======================================

This example opens an image stored in a FITS file and displays it to the screen.

This example uses `astropy.utils.data` to download the file, `astropy.io.fits` to open
the file, and `matplotlib.pyplot` to display the image.

-------------------

*By: Lia R. Corrales, Adrian Price-Whelan, Kelle Cruz*

*License: BSD*

-------------------

"""

##############################################################################
# Set up matplotlib and use a nicer set of plot parameters

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import glob
import os
SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd() + '/deblend0.05'
THRESH = 3

namecount = 0
def namegen():
    global namecount
    namecount += 1
    return DESTDIR + "\galaxyimage" + str(namecount) + ".png"

filenames = sorted(glob.glob(SOURCEDIR + '/psc*'))

import random
indices = []
for i in range(50):
    print
    indices.append(int(len(filenames)*random.random()))
    
for i in indices: #['C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc090275.3.fits']:#
    #print filename
    filename = filenames[i]
    print(str(namecount + 1) + ' ' + filename)
    #image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
    with fits.open(filename) as image_file:#'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc000014.6.fits')
        image_data = image_file[0].data
    

    #print(image_data.shape)
    
    '''
    plt.figure()
    plt.imshow(image_data, cmap='gray')
    plt.colorbar()
    '''
    #print("Done with first part")
    
    swappeddata = image_data.byteswap(inplace=True).newbyteorder()
    import sep
    bkg = sep.Background(swappeddata)
    
    data_sub = swappeddata# - bkg
    '''
    plt.figure()
    plt.imshow(data_sub, cmap='gray')#, extent=(0, numcols-0.5, numrows-0.5, -0.5))
    plt.colorbar()
    '''
    
    
    objects = sep.extract(data_sub, THRESH, err=bkg.globalrms, minarea = 5, deblend_cont = 0.05)
    from matplotlib.patches import Ellipse
    
    # plot background-subtracted image
    fig, ax = plt.subplots()
    import numpy as np
    m, s = np.mean(data_sub), np.std(data_sub)
    #im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
    #               vmin=m-s, vmax=m+s, origin='lower')
    
    im = ax.imshow(data_sub, interpolation='nearest', cmap='Greys_r', vmin = 0, vmax = 700)
    
    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
    plt.savefig(namegen(), dpi=150)
    
    im = ax.imshow(data_sub, interpolation='nearest', cmap='Greys_r')
    
    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
    plt.savefig(namegen(), dpi=150)
    
    '''
    flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
                                         3.0, err=bkg.globalrms, gain=1.0)
    
    for i in range(10):
        print("object {:d}: flux = {:f} +/- {:f}".format(i, flux[i], fluxerr[i]))
        
        
        add kron radius
        adjust magnitude by looking up stars magnitudes and making them match
        check for not in center, flag
        
    '''
