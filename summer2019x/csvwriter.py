# -*- coding: utf-8 -*-
"""
https://docs.python.org/3/library/csv.html
# http://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html
# http://docs.astropy.org/en/stable/io/fits/


find todos
remove background removal?
"""
import numpy as np
#from matplotlib.patches import Ellipse
import os
import glob
import csv
import sep
#import matplotlib.pyplot as plt
#from astropy.visualization import astropy_mpl_style
#plt.style.use(astropy_mpl_style)
#from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

SOURCEDIR = glob.glob('C:/Users/Faith/Desktop/noey2019summer/ps1hosts/psc1*')
DESTDIR = os.getcwd() + "\galaxiestest"
FILLER = [None]*7
THRESHOLD = 1.5
MINAREA = 5
#TODO check the following
MAG_ZERO # for calculating magnitude from flux


# make header
HEADER =['ID']
perImageHeaders = ['x', 'y', 'RA', 'DEC', 'KronMag', 'Angle', 'Ellipticity']
for i in range(3,7):
    for val in perImageHeaders:
        HEADER.append(val + '_' + str(i))


filenames = sorted(glob.glob(SOURCEDIR + '\\psc*.fits'))

#helper function to extract the data from a file
def getData(filename, threshold=THRESHOLD):
#TODO does this need to be recursive? if so move the following out, if not maybe
    #get rid of helper fxn
    image_file = fits.open(filename)
    image_data = image_file[0].data

    #https://www.aanda.org/index.php?option=com_article&access=doi&doi=10.1051/0004-6361/201015362&Itemid=129#FD14
    
    #fix byte order
    swappedData = image_data.byteswap(inplace=True).newbyteorder()
    
    # subtracting out background
    bkg = sep.Background(swappedData)
    swappedData = swappedData - bkg 
    
    objects = sep.extract(swappedData, threshold, 
                          err=bkg.globalrms, minarea = MINAREA)
    
    # if no sources detected, decrease threshold and try again. If threshold 
    # reaches threshhold min, send filler data
 #TODO if not objects:
    
    
    '''Find most likely host galaxy'''
    kronrad, _krflag = sep.kron_radius(swappedData, objects['x'], objects['y'], 
                                       objects['a'], objects['b'], 
                                       objects['theta'], 6.0)
    flux, _fluxerr, _flag = sep.sum_ellipse(swappedData, objects['x'], 
                                            objects['y'], objects['a'], 
                                            objects['b'], objects['theta'], 
                                            2.5*kronrad, subpix=1)
    magnitude = -2.5 * np.log10(flux) + MAG_ZERO
    
    #size is half light radius
    size, _flag = sep.flux_radius(swappedData, objects['x'], objects['y'], 
                      6.*objects['a'], 0.5, normflux=flux, subpix=5)
        
    separation = np.sqrt(((objects['x']-120)**2) + ((objects['y']-120)**2))
    
    # Observed number density of galaxies brighter than magnitude M (From Berger 2010)
    sigma = 10 ** (0.33 * (magnitude - 24) - 2.44) / (0.33 * np.log(10))

    # Effective radius
    R_effective = np.sqrt(np.abs(separation) + 4 * np.abs(size) ** 2)

    # Probability of chance coincidence
    chanceCoincidence = 1 - np.exp(-np.pi * R_effective ** 2 * sigma)
#TODO why so tiny chanceCoincidences    
    
    bestCandidate = np.argmin(chanceCoincidence)
    
    ellipticity = 1- (objects['b'][bestCandidate]/objects['a'][bestCandidate])
    

#TODO RA and DEC  
    centerRa = image_file.header['CRVAL1']
    centerDec = image_file.header['CRVAL2']

    return [objects['x']['bestCandidate'], objects['y'][bestCandidate],
            ra, dec, magnitude[bestCandidate],
            objects['theta'][bestCandidate], ellipticity]


# extract transient id# from filename
def getIdNum(filename):
    return filename[3:9]

def getFilter(filename):
    return filename[10:11]

def padAnyMissing(filename, lastIdNum, lastFilter):
#TODO
    #write 0s, also somehow write and reset if missing 3
    raise Exception("did not write padAnyMissing function")


with open(DESTDIR + "\galaxiesdata.csv", "w+") as destfile:
    csvwriter = csv.writer(destfile)
    csvwriter.writerow(HEADER)
    lastIdNum = None
    lastFilter = None
    cache = []
    for filename in filenames:
        idNum = getIdNum(filename)
        filterNum = getFilter(filename)
        # check if any missing files between current and previous, fill
        padAnyMissing(idNum, filterNum, lastIdNum, lastFilter)
        cache.append(getData(filename))
        
        # check if this is the last filter(6) of transient, if so write the 
        #cache to create line for this trasient. reset cache 
        if filterNum == '6':
            if len(cache) != 29:
                raise Exception("Incorrect row length")
            csvwriter.writerow(cache)
            cache = []
        
        lastIdNum = idNum
        lastFilter = filterNum
        
#TODO    missing padding

            
            
'''magnitude = -2.5 log(flux) + zeropoint

counts/exposuretime +25

480168.6 gives 18

to check, serach event, go to ned, closest event, magnitude in g band to see if close.
magnitude = -2.5 log(flux) + zeropoint

counts/exposuretime +25

480168.6 gives 18

to check, serach event, go to ned, closest event, magnitude in g band to see if close.
'''