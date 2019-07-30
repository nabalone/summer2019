# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:45:03 2019

@author: Faith
"""
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
DESTDIR = os.getcwd() + "\galaxiestest"
FILLER_VAL = None
THRESHOLD = 3
MINAREA = 5
MAG_ZERO =25 # for calculating magnitude from flux
DEBLEND_CONT = 0.05 # for sep.extract. 1.0 to turn of deblending (works best),
 # 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as

SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"

zs = []
ms = []

def getMag(filename):
    
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

    # objects to be ignored in host detection
    blacklist = set()
    
    '''Search SDSS for stars and blacklist them'''
#TODO remove
    realmags = []
#TODO fix
    for i in range(len(objects)):
#TODO change mindist to kron radius converted to degrees?
#TODO if too slow, change to make only one query
#TODO adjust mindist
        sdssTable = SDSS.query_region(coordinates=objCoords[i], radius=MINDIST, 
                              photoobj_fields=['ra','dec','type', 'mode', 'i'])
        ''' IMPORTANT fields are here:
            http://skyserver.sdss.org/dr12/en/help/browser/browser.aspx#&&history=description+PhotoObjAll+U
        '''

        #if table of objects at that location exists and 
        # a star and no galaxies are there. star is type 6, galaxy is type 3
        #if sdssTable and 6 in sdssTable['type'] and not 3 in sdssTable['type']:
        # nevermind, actually if the primary object (mode=1) is star then blacklist
        realmag = None
        if sdssTable:    
            for j in range(len(sdssTable)):
                if sdssTable['mode'][j] == 1 and sdssTable['type'][j] == 6:
                    blacklist.add(i)
#TODO remove
                    realmag = sdssTable['i'][j]
                    break
        realmags.append(realmag)
        # remove any 'nan' kronrads and blacklist the objects
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
        # rough convert to arcsecs
        # TODO refine?
        degPerPix = image_file[0].header['CDELT1']
        size = size * degPerPix *3600
          
        separation = objCoords.separation(eventCoords).arcsec
        #np.sqrt(((objects['x']-120)**2) + ((objects['y']-120)**2))
        
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


