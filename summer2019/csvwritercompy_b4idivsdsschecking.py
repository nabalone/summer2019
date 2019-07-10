# -*- coding: utf-8 -*-

"""
find todos
remove background removal?
"""
import numpy as np
import os
import glob
import random
import csv
import sep
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy.wcs import WCS

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.sdss import SDSS

SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd() + "\galaxiestest"
FILLER_VAL = None
THRESHOLD = 3
MINAREA = 5
MAG_ZERO =25 # for calculating magnitude from flux
DEBLEND_CONT = 0.05 # for sep.extract. 1.0 to turn of deblending (works best),
 # 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as

#For debugging:
WRITE_CSV = "\galaxiesdata.csv" # filename to write to or None
PRINT_DATA = False
CHECK_DISTANCE = None #print all files with most likely host farther than this int
PLOT_ALL = True
PLOT_ERR = False #plots all that give errors or low probability
PLOT_DIR = os.getcwd() + '/allfounds2' # where to put plot images
ONLY_FLAG_ERRORS = True # catch errors, print filename, move on
FILES = 'specified' #options are 'all', 'preset random', 'new random', 'range', 'specified'
SPECIFIED = [SOURCEDIR + '\psc450082.3.fits']
RANGE = (200, 210)

# make header
HEADER =['ID']
perImageHeaders = ['KronRad', 'separation', 'x', 'y', 'RA', 'DEC', 'KronMag', 'Angle', 'Ellipticity']
for i in range(3,7):
    for val in perImageHeaders:
        HEADER.append(val + '_' + str(i))

#figure out which files to use based on value specified at top
if FILES == 'all' or FILES =='range' or FILES =='new random':
    filenames = sorted(glob.glob(SOURCEDIR + '/psc*.[3-6].fits'))
    if FILES == 'range':
        filenames = filenames[RANGE[0]:RANGE[1]]
elif FILES == 'preset random':
    from presetrandomfiles import fileset
    filenames = fileset
elif FILES == 'new random':
    temp = []
    for i in range(100):
        temp.append(filenames[int(len(filenames)*random.random())])
    filenames = temp
elif FILES == 'specified':
    filenames = SPECIFIED
else:
    raise Exception('invalid FILE specification')
    
# for naming plots files
namecount = 0
def namegen():
    global namecount
    namecount += 1
    return PLOT_DIR + "\galaxyimage" + str(namecount) + ".png"
    
def plot(data, objects, blacklist, bestCandidate, myVmin=None, myVmax=None):
    # make the destination directory if it does not exist
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)
    to_circle = range(len(objects))
    
    fig, ax = plt.subplots()

    if myVmax == None or myVmin == None:    
        _im = ax.imshow(data, interpolation='nearest', cmap='gray')
    else:
        _im = ax.imshow(data, interpolation='nearest', cmap='gray',
                        vmin = myVmin, vmax = myVmax)
        
    # plot an ellipse for each object
    for i in to_circle:
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        if i==bestCandidate:
            e.set_edgecolor('green')
        elif i in blacklist:
            e.set_edgecolor('blue')
        else:
            e.set_edgecolor('red')
        ax.add_artist(e) 
    plt.savefig(namegen(), dpi=150)    



'''main'''
def errorProtocol():
    print("error: " + filename)
    if PLOT_ERR:
        plot(swappedData, objects, bestCandidate, myVmin = 0, myVmax = 1000)
        plot(swappedData, objects, bestCandidate)
    if ONLY_FLAG_ERRORS:
        return
    else:
        raise
    
if WRITE_CSV:
    #create destination directory if it does not exist
    if not os.path.isdir(DESTDIR):
        os.mkdir(DESTDIR)
    destfile = open(DESTDIR + WRITE_CSV, "w+")
    csvwriter = csv.writer(destfile)
    csvwriter.writerow(HEADER)
if PRINT_DATA:
    print HEADER

# idnum of first file. See below for extraction logic
lastIdNum = int(filenames[0].split('.')[-2])
lastFilter = 2
cache = []

for filename in filenames:
    print filename
    # to extract transient and filter number from filename of the form
    # */psc190369.6.fits
    dotSplit = filename.split('.')
    # transient identifier is everything after the last 'c' but before 
    # the second to last dot
    idNum = int(dotSplit[-3].split('c')[-1]) 
    # filter number is between second to last dot and last dot
    filterNum = int(dotSplit[-2]) # filter of this specific image
    
    # if you are on a new event numbeer, 
    # fill the old row, write it, reset filter to 2, reset cache
    if lastIdNum != idNum:
        cache.extend([FILLER_VAL]*(29-len(cache)))
        if WRITE_CSV:
            csvwriter.writerow(cache)
        if PRINT_DATA:
            print cache
        cache = [idNum]
        lastFilter = 2
        
    # fill any missing data in current row
    for i in range(filterNum - lastFilter - 1):
        cache.extend([FILLER_VAL]*len(perImageHeaders))

    '''start getData block'''
#TODO does this need to be recursive on trying threshholds?
    image_file = fits.open(filename)
    image_data = image_file[0].data
    
    #fix byte order
    swappedData = image_data.byteswap(inplace=True).newbyteorder()
    
    # subtracting out background
    bkg = sep.Background(swappedData)
    if SUBTRACT_BACKGROUND:
        swappedData = swappedData - bkg 
    
    objects = sep.extract(swappedData, THRESHOLD, err=bkg.globalrms, 
                          minarea = MINAREA, deblend_cont = DEBLEND_CONT)        
#TODO if not objects:

    '''Find most likely host galaxy'''
    # how to calculate kron radius and flux from 
    # https://sep.readthedocs.io/en/v1.0.x/apertures.html
    kronrad, _krflag = sep.kron_radius(swappedData, objects['x'], objects['y'], 
                                       objects['a'], objects['b'], 
                                       objects['theta'], 6.0)
    
    w = WCS(filename)
    # the last argument is 1 b/c first pixel of FITS should be 1,1 not 0,0?
    ra, dec = w.all_pix2world(objects['x'], objects['y'], 1)
    centerRa, centerDec = w.all_pix2world(120, 120, 1)
    # get rid of the array encasing
    centerRa = centerRa.item(0)
    centerDec = centerDec.item(0)

    '''Search SDSS for stars and blacklist them'''
    center = SkyCoord(centerRa*u.deg, centerDec*u.deg, frame='icrs')
    
    # to get image dimensions in wcs:
    maxX = image_file[0].header['NAXIS1']
    maxY = image_file[0].header['NAXIS2']
    ra1, dec1 = w.all_pix2world(1,1,1)
    raMaxX, decMaxY = w.all_pix2world(maxX, maxY, 1)
    imageWidth = ra1 - raMaxX
    imageHeight = dec1 - decMaxY

    #WHY WON'T IT LET ME SPECIFY WIDTH???
    # use circle circumscribing image until I can figure out how to use width
    rad = (np.sqrt(2) * imageWidth / 2.) * u.deg
    sdssTable = SDSS.query_region(coordinates=center, radius=rad, 
                                  photoobj_fields=['ra','dec','type', 'mode'])
                                     #width=imageWidth*u.deg,height=imageHeight)
    objCoords = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs') 
    
    # objects to be ignored in host detection
    blacklist = set()
    try:
        sdssCoords = SkyCoord(sdssTable['ra']*u.deg, sdssTable['dec']*u.deg)
        indices,dists,_dists = objCoords.match_to_catalog_sky(sdssCoords)
        # TODO check if this all works. no sdss objs on 470047.3 except one off field.

        
        
        for i in range(len(indices)):
            # object i is matched to a star in sdss within min dist
            if dists[i] < MINDIST and sdssTable[indices[i]]['type'] == 6: #obj is a star
                blacklist.add(i)
    except TypeError: # no matching stars found
        pass

    # remove any 'nan' kronrads and blacklist the objects
    for i in range(len(objects)):
        #check if ki is nan:
        if (not kronrad[i] == 0) and (not kronrad[i] < 0) and (not kronrad[i] > 0):
            # remove that object's column
            kronrad[i] = 0.1
            blacklist.add(i)
    try:
        flux, _fluxerr, _flag = sep.sum_ellipse(swappedData, objects['x'], 
                                                objects['y'], objects['a'], 
                                                objects['b'], objects['theta'], 
                                                2.5*kronrad, subpix=1)


        magnitude = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + MAG_ZERO 
        
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
#TODO check reasonable chance coincidence 
#TODO exclude stars from catalog
        # exclude any blacklisted
        for i in blacklist:
            chanceCoincidence[i] = 1
        bestCandidate = np.argmin(chanceCoincidence)
        
        if CHECK_DISTANCE:
            if separation[bestCandidate] > CHECK_DISTANCE:
                print("far: " + filename + " "+ str(separation[bestCandidate]))
        ellipticity = 1 - (objects['b'][bestCandidate]/objects['a'][bestCandidate])

        image_file.close()
            
        cache.extend([kronrad[bestCandidate], separation[bestCandidate], 
                objects['x'][bestCandidate], objects['y'][bestCandidate],
                ra, dec, magnitude[bestCandidate],
                objects['theta'][bestCandidate], ellipticity])
        ''' end getData block'''            
        
        lastIdNum = idNum
        lastFilter = filterNum
        
        if PLOT_ALL:
            plot(swappedData, objects, blacklist, bestCandidate)
            plot(swappedData, objects, blacklist, bestCandidate, 
                 myVmin = 0, myVmax = 1000)
            
    except:
        errorProtocol()
    
# write last data
cache.extend([FILLER_VAL]*(29-len(cache)))
if WRITE_CSV:
    csvwriter.writerow(cache)
    destfile.close()
if PRINT_DATA:
    print cache
        
        
# =============================================================================
# if __name__ == "__main__":
#     main()
# 
# =============================================================================
