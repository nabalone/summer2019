# -*- coding: utf-8 -*-

"""
find todos
remove background removal?
"""
import numpy as np
import pandas as pd
import os
import glob
import random
import csv
import sep
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from matplotlib.patches import Ellipse, RegularPolygon
from astropy.io import fits
from astropy.wcs import WCS

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.sdss import SDSS

SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd() + "\magnitudes"
FILLER_VAL = None
THRESHOLD = 3
MINAREA = 5
MAG_ZERO =25 # for calculating magnitude from flux
DEBLEND_CONT = 0.05 # for sep.extract. 1.0 to turn of deblending (works best),
 # 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as

#For debugging:
WRITE_CSV = "\magdata.csv" # filename to write to or None

ONLY_FLAG_ERRORS = False # catch errors, print filename, move on
FILES = 'new random' #options are 'all', 'preset random', 'new random', 'range', 'specified', 'nonsquare
SPECIFIED = [SOURCEDIR + '\psc050601.3.fits']
#[SOURCEDIR + '\psc000010.3.fits', SOURCEDIR + '\psc000010.4.fits', SOURCEDIR + '\psc000010.5.fits', SOURCEDIR + '\psc000010.6.fits']
RANGE = (200, 210)

# make header
HEADER = ['ID', 'Filter', 'obj', 'mymag', 'realmag', 'dif']#z', 'i', 'r', 'g']
#figure out which files to use based on value specified at top
if FILES == 'all' or FILES =='range' or FILES =='new random':
    filenames = sorted(glob.glob(SOURCEDIR + '/psc*.[3-6].fits'))
    if FILES == 'range':
        filenames = filenames[RANGE[0]:RANGE[1]]
    if FILES == 'new random':
        temp = []
        for i in range(10):
            temp.append(filenames[int(len(filenames)*random.random())])
        filenames = temp
elif FILES == 'preset random':
    from presetrandomfiles import fileset
    filenames = fileset[:12]
elif FILES == 'nonsquare':
    from nonsquare import fileset
    filenames = fileset
elif FILES == 'specified':
    filenames = SPECIFIED
else:
    raise Exception('invalid FILE specification')

 
if not os.path.isdir(DESTDIR):
    os.mkdir(DESTDIR)
destfile = open(DESTDIR + WRITE_CSV, "w+")
csvwriter = csv.writer(destfile)
csvwriter.writerow(HEADER)

cache = []

for filename in filenames:
    '''
    for filternum in range(3,7):
        filename = origfilename[:-6] + str(filternum) + origfilename[-5:]      
    '''
    '''start getData block'''
    
        # to extract transient and filter number from filename of the form
    # */psc190369.6.fits
    dotSplit = filename.split('.')
    # transient identifier is everything after the last 'c' but before 
    # the second to last dot
    idNum = int(dotSplit[-3].split('c')[-1]) 
    # filter number is between second to last dot and last dot
    filterNum = int(dotSplit[-2]) # filter of this specific image

    image_file = fits.open(filename)
    image_data = image_file[0].data
    
    #fix byte order
    swappedData = image_data.byteswap(inplace=True).newbyteorder()
    
    # subtracting out background
    bkg = sep.Background(swappedData)
    if SUBTRACT_BACKGROUND:
        swappedData = swappedData - bkg 
        
    # get noise data
    '''
    noisefilename = filename[:-4] + 'noise.' + filename[-4:]
    with fits.open(noisefilename) as noise_file:
        noise_data = noise_file[0].data
    noise_data2 = noise_data.byteswap(inplace=True).newbyteorder()
    '''

    objects = sep.extract(swappedData, THRESHOLD,  err=bkg.globalrms, #var=noise_data2 ,
                          minarea = MINAREA, deblend_cont = DEBLEND_CONT)        
#TODO if not objects:

    '''Find most likely host galaxy'''
    # how to calculate kron radius and flux from 
    # https://sep.readthedocs.io/en/v1.0.x/apertures.html
    kronrad, _krflag = sep.kron_radius(swappedData, objects['x'], objects['y'], 
                                       objects['a'], objects['b'], 
                                       objects['theta'], 6.0)
    for i in range(len(objects)):
        # remove any 'nan' kronrads and blacklist the objects
        #check if ki is nan:
        if (not kronrad[i] == 0) and (not kronrad[i] < 0) and (not kronrad[i] > 0):
            # remove that object's column
            kronrad[i] = 0.001

    flux, _fluxerr, _flag = sep.sum_ellipse(swappedData, objects['x'], 
                                            objects['y'], objects['a'], 
                                            objects['b'], objects['theta'], 
                                            2.5*kronrad, subpix=1)


    magnitude = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + MAG_ZERO 

    w = WCS(filename)
    # the last argument is 1 b/c first pixel of FITS should be 1,1 not 0,0?
    ra, dec = w.all_pix2world(objects['x'], objects['y'], 1)
    objCoords = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs') # coords of all objs
    
    # find pixel coordinates of event
    # TODO check if this is unnecessary for squares
    db = pd.read_table('alertstable_v3',sep=None,index_col = False, 
                       engine='python')
    db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False, 
                        engine='python')
    db = db.append(db2,ignore_index=True)
    event = db.where(db['eventID'] == idNum).dropna()
    eventRa = event['ra'].values[0] #values gives np arrays
    eventDec = event['dec'].values[0] #'hh:mm:ss.sss'
    
    # converting to degrees
    eventCoords = SkyCoord(eventRa, eventDec, unit=(u.hourangle, u.deg))
    eventRa = eventCoords.ra.deg
    eventDec = eventCoords.dec.deg
#TODO: ircs or fk5
    
    # get event pixel coords
    eventX, eventY = w.all_world2pix(eventRa, eventDec, 1)
    
    '''Search SDSS for stars and blacklist them'''
#TODO remove
    realmags = []
#TODO fix
    for i in range(len(objects)):
#TODO change mindist to kron radius converted to degrees?
#TODO if too slow, change to make only one query
#TODO adjust mindist
        sdssTable = SDSS.query_region(coordinates=objCoords[i], radius=MINDIST, 
                              photoobj_fields=['ra','dec','type', 'mode', 
                                               'modelMag_z', 'modelMag_i', 'modelMag_r', 'modelMag_g'])
        ''' IMPORTANT fields are here:
            http://skyserver.sdss.org/dr12/en/help/browser/browser.aspx#&&history=description+PhotoObjAll+U
        '''

        #if table of objects at that location exists and 
        # a star and no galaxies are there. star is type 6, galaxy is type 3
        #if sdssTable and 6 in sdssTable['type'] and not 3 in sdssTable['type']:
        # nevermind, actually if the primary object (mode=1) is star then blacklist
        if sdssTable:    
            for j in range(len(sdssTable)):
                if sdssTable['mode'][j] == 1:
                    cache.append(idNum)
                    cache.append(filterNum)
                    cache.append(i)
                    cache.append(magnitude[i])
                    if filterNum == 3:
                        realmag = sdssTable[j]['modelMag_z']
                    elif filterNum == 4:
                        realmag = sdssTable[j]['modelMag_i']
                    elif filterNum == 5:
                        realmag = sdssTable[j]['modelMag_r']
                    elif filterNum == 6:
                        realmag = sdssTable[j]['modelMag_g']
                    
                    cache.append(realmag)
                    print filterNum
                    print realmag

                    cache.append(magnitude[i] - float(realmag))
                    
#                    cache.append(sdssTable[j]['modelMag_z'])
#                    cache.append(sdssTable[j]['modelMag_i'])
#                    cache.append(sdssTable[j]['modelMag_r'])
#                    cache.append(sdssTable[j]['modelMag_g'])
                    csvwriter.writerow(cache)
                    cache = []
                    break

    image_file.close()

destfile.close()

        
        
# =============================================================================
# if __name__ == "__main__":
#     main()
# 
# =============================================================================
