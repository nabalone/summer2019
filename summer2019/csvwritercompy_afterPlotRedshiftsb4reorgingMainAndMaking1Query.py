# -*- coding: utf-8 -*-

#TODO: all constants at top
#TODO: vertical line
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
WRITE_CSV = None #"\galaxiesdata.csv" # filename to write to or None
PRINT_DATA = False
CHECK_DISTANCE = None #print all files with most likely host farther than this int
PLOT_ALL = False
PLOT_ERR = False #plots only files that give errors or low probability
PLOT_DIR = os.getcwd() + '/allfounds' # where to put plot images
ONLY_FLAG_ERRORS = False # catch errors, print filename, move on
FILES = 'specified' #options are 'all', 'preset random', 'new random', 'range', 'specified', 'nonsquare
SPECIFIED = [SOURCEDIR + '\psc050601.4.fits']
#[SOURCEDIR + '\psc000010.3.fits', SOURCEDIR + '\psc000010.4.fits', SOURCEDIR + '\psc000010.5.fits', SOURCEDIR + '\psc000010.6.fits']
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
    if FILES == 'new random':
        temp = []
        for i in range(10):
            temp.append(filenames[int(len(filenames)*random.random())])
        filenames = temp
elif FILES == 'preset random':
    from presetrandomfiles import fileset
    getData(fileset[:12])
elif FILES == 'nonsquare':
    from nonsquare import fileset
    filenames = fileset
elif FILES == 'specified':
    filenames = SPECIFIED
else:
    raise Exception('invalid FILE specification')
    
''' My stupid workaround for debugging, so I can access the variables from the console'''
objects=None
bestCandidate=None
magnitudes=None
sdssTable=None
chanceCoincidence=None
separation=None
filename=None
blacklist=None
    
    
# for naming plots files
namecount = 0
def namegen():
    global namecount
    namecount += 1
    return PLOT_DIR + "\galaxyimage" + str(namecount) + ".png"
    
def plot(data, objects, blacklist, bestCandidate, myVmin=None, myVmax=None,
         myEventX = None, myEventY = None):
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
    # triangle on event location
    p = RegularPolygon((int(myEventX), int(myEventY)), 3, radius=3)
    p.set_edgecolor('purple')
    p.set_facecolor('purple')
    ax.add_artist(p)
        
    # plot an ellipse for each object
    for i in to_circle:
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
# TODO check if above angle conversion is correct. Where from?
        e.set_facecolor('none')
        if i==bestCandidate:
            e.set_edgecolor('green')
        elif i in blacklist:
            e.set_edgecolor('blue')
        else:
            e.set_edgecolor('red')
        ax.add_artist(e) 
    plt.savefig(namegen(), dpi=150)    
    
zs = []
ms = []
def plotRedshifts():
    zfile = open('new_ps1z.dat', 'r')
    zfile.readline() #get rid of header
    for i in range(50):
        line = zfile.readline()
        parts = line[:-1].split('\t') #strip line ending and split by tabs
        eventId = parts[0][3:]
        redshift = float(parts[1])
        filename = SOURCEDIR + '/psc' + eventId + '.6.fits'
        try:
            mag = main([filename], special='mag_only')
            ms.append(mag)
            zs.append(redshift)
        except IOError as e:
            print(e)
    zfile.close()
    
    plt.plot(zs, ms, 'bo')
    plt.savefig('redshift_vs_magnitude', dpi=150)
 
def main(filenames, special=None):
    ''' My stupid workaround for debugging, so I can access the variables from the console'''
    global objects
    global bestCandidate
    global magnitudes
    global sdssTable
    global chanceCoincidence
    global separation
    global filename
    global blacklist
    
    def errorProtocol():
        print("error: " + filename)
        if PLOT_ERR:
            plot(swappedData, objects, bestCandidate, myVmin = 0, myVmax = 1000, 
                 myEventX = eventX, myEventY = eventY)
            plot(swappedData, objects, bestCandidate, 
                 myEventX = eventX, myEventY = eventY)
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
                
            cache.extend([kronrad[bestCandidate], separation[bestCandidate], 
                    objects['x'][bestCandidate] - eventX, objects['y'][bestCandidate] - eventY,
                    ra[bestCandidate], dec[bestCandidate], magnitude[bestCandidate],
                    objects['theta'][bestCandidate], ellipticity])
            ''' end getData block'''            
            
            lastIdNum = idNum
            lastFilter = filterNum
            
    #        print filename
    #        print objCoords.to_string('hmsdms', sep=':')
    #        print magnitude
    #        print ''
    
            if special=='mag_only':
                return magnitude[bestCandidate]
            
            if PLOT_ALL:
                plot(swappedData, objects, blacklist, bestCandidate,
                     myEventX = eventX, myEventY = eventY)
                plot(swappedData, objects, blacklist, bestCandidate, 
                     myVmin = 0, myVmax = 1000, myEventX = eventX, myEventY = eventY)
                
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
