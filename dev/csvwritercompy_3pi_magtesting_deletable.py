# -*- coding: utf-8 -*-

#TODO: all constants at top
#TODO: vertical line
'''
fits files must be named pscxxxxxx.f.fits
and FILTERS must be such that filter number f corresponds with filter filters[f]
alertstable_vs and alertstablevs.lasthalf must be in directory for sn locations
sdssTableGroupIndex.py, sdss_queries.dat, new_ps1z.dat, 'ps1confirmed_only_sne_without_outlier.txt',
hosts.dat

'''

import numpy as np
import pandas as pd
import os
import glob
import random
import csv
import sep
import ast
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from matplotlib.patches import Ellipse, RegularPolygon
from astropy.io import fits
from astropy.wcs import WCS

from astropy import units as u
from astropy.coordinates import SkyCoord
#from astroquery.sdss import SDSS
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table#, vstack
import json
#from sdssTableGroupIndex import sdssTableGroupIndex

ERRORFILE = 'errorfile2.txt'
SOURCEDIR = os.getcwd()
DESTDIR = os.getcwd()
FILLER_VAL = None
THRESHOLD = 3.
MINAREA = 4
DEBLEND_CONT = 0.01 # for sep.extract. 1.0 to turn off deblending, 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as
FILTERS = [None, None, None, 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z']
LIKELIHOOD_THRESH = 0.2
TYPES = ['SNIIn', 'SNIa', 'SNII', 'SNIbc', 'SLSNe']
TYPE_COLORS = {'SNIIn':'co', 'SNIa':'ro', 'SNII': 'bo', 'SNIbc':'go', 'SLSNe': 'mo'}
LOWEST_MAG = 40
# meds used for m_0 if no sdss stars in image, or if m_0 is an outlier
# meaning outside of upper/lower bounds
#FILTER_M0s = [None, None, None,
#              {'upper': 23.625, 'default': 23.125, 'lower': 22.625}, # filter 3
#              {'upper': 23.4, 'default': 22.9, 'lower': 22.4}, # filter 4
#              {'upper': 23.4, 'default': 22.9, 'lower': 22.4}, # filter 5
#              {'upper': 24.125, 'default': 23.625, 'lower': 23.125}] # filter 6

#update from first 120 files:
FILTER_M0s = [None, None, None,
              {'upper': 23.41, 'default': 22.918, 'lower': 22.41}, # filter 3
              {'upper': 23.32, 'default': 22.822, 'lower': 22.32}, # filter 4
              {'upper': 24.15, 'default': 23.652, 'lower': 23.15}, # filter 5
              {'upper': 24.04, 'default': 23.540, 'lower': 23.04}] # filter 6


#USAGE FLAGS:
USE_3PI = True
WRITE_CSV = "\galaxiesdata2.csv" # filename to write to or None
MAG_TEST_ALL = False
MAG_TEST_STDEV = False
PLOT_REDSHIFTS = False

MAGFILE = 'newmagtest.csv'
#Do not write csv an mag test simultaneously because there is only 1 writer
#TODO separate the csvwriters
#ACTUALLY< REMOVING MAG WRITER
#if MAG_TEST_ALL or MAG_TEST_STDEV:
#    WRITE_CSV = False
PRINT_DATA = True

CHECK_DISTANCE = 5 #print all files with most likely host farther than this arcsecs
PLOT_ALL = False
PLOT_ERR =  False #plots only files that give errors or low probability
PLOT_DIR = os.getcwd() + '/plots2' # where to put plot images
ONLY_FLAG_ERRORS = False # catch errors, print filename, move on
FILES = 'range' #options are 'all', 'preset random', 'new random', 'range', 'specified', 'nonsquare

#TODO delete
SPECIFIED = []
to_check = [160103, 180313, 590123, 50296, 90034, 50601]
for f in to_check:
    SPECIFIED.extend(glob.glob((SOURCEDIR + '/ps1hosts/psc*%i*.[3-6].fits' % f)))

SPECIFIED = [SOURCEDIR + '/3pi_hosts/ps_000006.3.fits']
RANGE = (0, 100)
m0collector = [None, None, None, [], [], [], []]

'''make header'''
HEADER =['ID']
#perImageHeaders = ['KronRad', 'separation', 'x', 'y', 'RA', 'DEC', 'KronMag', 'Angle', 'Ellipticity']
perImageHeaders = ['KronRad (kpc)', 'separation (kpc)', 'area (kpc^2)', 'sep/sqrt(area) (kpc)',
                   'x', 'y','KronMag', 'Abs. Mag', 'Angle',
                   'Ellipticity', 'RA', 'Host RA', 'DEC', 'Host Dec',
                   'Discrepency (arcsecs)', 'Z', 'SDSS Photoz', 'pixelRank', 'chance coincidence']
for i in range(3,7):
    for val in perImageHeaders:
        HEADER.append(val + '_' + str(i))

# remove error file
if os.path.exists(ERRORFILE):
  os.remove(ERRORFILE)

# for naming plots files
namecount = 0
def namegen():
    global namecount
    namecount += 1
    return PLOT_DIR + "\galaxyimage" + str(namecount) + ".png"

def plot(data, objects, blacklist, bestCandidate, myVmin=None, myVmax=None,
         myEventX = None, myEventY = None, special_plot_location=None, title=''):
    # make the destination directory if it does not exist
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)
        
    special_plot_dir = ''
    if special_plot_location:
        special_plot_dir = PLOT_DIR + special_plot_location + '.png'
            
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
    plt.title(title)
    genname = namegen()
    name = special_plot_dir if special_plot_dir else genname
    plt.savefig(name, dpi=150)
    plt.show()
    plt.close()

def getPixelRank(swappedData, eventX, eventY, segmap, bestCandidate):
    global a
    global pixels
    global location
    a = np.where(segmap == bestCandidate+1)
    pixels = []
    for k in range(len(a[0])):
        pixels.append((a[0][k], a[1][k])) #list of tuples of coords
    #pixels.append((int(eventX), int(eventY))) #in case event is outside host object
    def sortkey(x):
        a, b = x
        return swappedData[a][b]
    #if pixel not in object
    if not (int(eventX), int(eventY)) in pixels:
        return 0
    pixels.sort(key = sortkey)
    location = pixels.index((int(eventX), int(eventY)))
    return float(location)/float(len(pixels))

def FWHM(arr_y):
    difference = max(arr_y) - min(arr_y)
    HM = difference / 2
    
    peak = arr_y.argmax()  # or in your case: arr_y.argmin()
    
    left = (np.abs(arr_y[:peak] - HM)).argmin()
    right = peak + (np.abs(arr_y[peak:] - HM)).argmin()

    return right - left


def extraction(filenames, additional=False):
    ''' My stupid workaround for debugging, so I can access the variables from the console'''
    global objects
    global bestCandidate
    global magnitude
    global sdssTable
    global chanceCoincidence
    global separation
    global filename
    global blacklist
    global swappedData
    global image_file
    global ra
    global dec
    global hostRa
    global hostDec
    global eventRa
    global eventDec
    global colRealMags
    global colMyMags
    global colFluxes
    global zdict
    global all_myMags
    global all_realMags
    global BAD_IMAGES
    global GOOD_IMAGE
    global m_0s
    global m_0
    global size
    global R_effective
    global hostsData
    global all_myMagsFiltered
    global all_realMagsFiltered
    global all_redshifts
    global all_kronMags
    global all_kronRads
    global kronrad
    global flux
    global z
    global photozs
    global eventX
    global eventY
    global BAD_IMAGES
    global GOOD_LOCATIONS
    global GOOD_PHOTOZS
    global BAD_TO_GOOD
    global GOOD_PHOTOZS_TO_BAD
    global cache


    all_realMags = [[], [], [], [], [], [], []]
    all_myMags = [[], [], [], [], [], [], []]
    GOOD_IMAGE = None
    BAD_IMAGES = []

    all_redshifts = {}
    all_kronMags = {}
    all_kronRads = {}
    for snType in TYPES:
        all_redshifts[snType] = []
        all_kronMags[snType] = []
        all_kronRads[snType] = []
    global errorProtocol
    def errorProtocol(e, specified_file=None, green=None):
        if specified_file:
            curFile = specified_file
        else:
            curFile = filename
#TODO fix
#        try:
#            errorString = str(namecount) + '\n' + str(e) + " " + curFile + " dist: " + \
#                          str(separation[bestCandidate]) + " chanceCoincidence: " + \
#                          str(chanceCoincidence[bestCandidate])  + '\n'
#        except: #probably separation or chanceCoincidence not calculated
#            print("x")
        errorString = str(namecount) + '\n' + str(e) + " " + curFile + '\n'
        print(errorString)
        with open(ERRORFILE, 'a+') as errorfile:
            errorfile.write(errorString)

        chosen = green if green else bestCandidate

        special_plot_location = '/NoGood/' + e if e[:5] == "Best." else None
        special_plot_location2 = '/NoGood/' + e + 'b' if e[:5] == "Best." else None
        if not os.path.isdir(PLOT_DIR + '/NoGood/'):
            os.mkdir(PLOT_DIR + '/NoGood/')
        e = str(e)
        if PLOT_ERR:
            padded_e = e + '                                                                         '
            my_title = padded_e[:30] + '\n' + curFile[-16:]
            plot(swappedData, objects, blacklist, chosen, myVmin = 0, myVmax = 1000,
                 myEventX = eventX, myEventY = eventY, 
                 special_plot_location=special_plot_location, title=my_title)
            plot(swappedData, objects, blacklist, chosen,
                 myEventX = eventX, myEventY = eventY, 
                 special_plot_location=special_plot_location2, title=my_title)
#        if ONLY_FLAG_ERRORS or e=='far' or e=='unlikely':
#            return
#        else:
#            print("raising")
#            raise Exception("error caught")

    if additional:
        '''load data for the only in 3pi additions'''
        db_3pi = pd.read_csv('3pi/all_additions.csv')
        
        with open("sdss_queries_index_3pi.txt") as s:
            sdssTableGroupIndex = json.load(s)
         
        '''load prerun sdss queries'''
        fullSdssTable = Table.read('sdss_queries_3pi.dat', format='ascii')
        fullSdssTable = fullSdssTable.group_by('idnum')
    
    else:
        ''' load event redshifts' dictionary '''
        zdict = {}
        zfile = open('new_ps1z.dat', 'r')
        zfile.readline() #get rid of header

        for line in zfile:
            parts = line.split()
            eventId = parts[0][3:]
            redshift = float(parts[1])
            zdict[eventId] = redshift
        zfile.close()
        zdict['100014'] = 0.357
        zdict['300220'] = 0.094
        zdict['380108'] = 0.159
        

        '''load event type dictionary'''
        typeDict = {}
        typefile = open('ps1confirmed_only_sne_without_outlier.txt', 'r')
        typefile.readline() #get rid of header
        for line in typefile:
            parts = line.split()
            eventId = parts[0][3:]
            eventType = parts[1]
            typeDict[eventId] = eventType
        typefile.close()

        '''load event location dictionary'''
        db = pd.read_table('alertstable_v3',sep=None,index_col = False,
                       engine='python')
        db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False,
                        engine='python')
        db = db.append(db2,ignore_index=True)

        
        with open("sdss_queries_index.txt") as s:
            sdssTableGroupIndex = json.load(s)

        '''load prerun sdss queries'''
        fullSdssTable = Table.read('sdss_queries.dat', format='ascii')
        fullSdssTable = fullSdssTable.group_by('idnum')

#TODO check appending works ok
    if WRITE_CSV: #additionals always come after originals, so we can
    #write fresh incl. header if originals and append if additionals
        #create destination directory if it does not exist
        if not os.path.isdir(DESTDIR):
            os.mkdir(DESTDIR)
        if not additional:
            destfile = open(DESTDIR + WRITE_CSV, "w+")
            csvwriter = csv.writer(destfile)
            csvwriter.writerow(HEADER)
        else:
            destfile = open(DESTDIR + WRITE_CSV, "a+")
            csvwriter = csv.writer(destfile)            
    if PRINT_DATA:
        print(HEADER)

    # idnum of first file. See below for extraction logic
    lastIdNum = int(filenames[0].split('.')[-2])
    lastFilter = 2
    cache = []

#    if MAG_TEST_ALL or MAG_TEST_STDEV:
#        magdestfile = open(MAGFILE, "w+")
#        csvwriter = csv.writer(magdestfile)

    for filename in filenames:
        ''' clean up and set up '''
        #clear bestCandidate to avoid bad plotting
        bestCandidate = None

        # to extract transient and filter number from filename of the form
        # */psc190369.6.fits
        dotSplit = filename.split('.')
        # transient identifier is everything after the last 'c' but before
        # the second to last dot
#TODO check
        idNumString = dotSplit[-3][-6:]#.split('c')[-1] #used for getting hectospec data
        idNum  = int(idNumString)
        # filter number is between second to last dot and last dot
        filterNum = int(dotSplit[-2]) # filter of this specific image

        
        # if you are on a new event number,
        # fill the old row, write it, reset filter to 2, reset cache
        if lastIdNum != idNum:
            cache = [idNum]
            BAD_IMAGES = []
            GOOD_LOCATIONS = []
            GOOD_PHOTOZS = []
            BAD_TO_GOOD = []
            GOOD_PHOTOZS_TO_BAD = []
        '''
            for imageDict in BAD_IMAGES:
                # these are for an event for which we found no host in center:
                errorProtocol('far', specified_file=str(lastIdNum) + '.' + \
                              str(imageDict['filterNum']) + ' no error info',
                              green = imageDict['old_best'])
            GOOD_IMAGE = None
            BAD_IMAGES = []

            cache.extend([FILLER_VAL]*(1+4*len(perImageHeaders)-len(cache)))
            if WRITE_CSV:
                csvwriter.writerow(cache)
            if PRINT_DATA:
                print(cache)
            cache = [idNum]
            lastFilter = 2

        # fill any missing data in current row
        for i in range(filterNum - lastFilter - 1):
            cache.extend([FILLER_VAL]*len(perImageHeaders))
        '''
        image_file = fits.open(filename)


        if not additional:
            ''' get data on the image '''
            # find pixel coordinates of event
            event = db.where(db['eventID'] == idNum).dropna()
            eventRa = event['ra'].values[0] #values gives np arrays
            eventDec = event['dec'].values[0] #'hh:mm:ss.sss'
        
        else:
            event = db_3pi.iloc[int(idNum)]
            #event = db_3pi.where(db_3pi[list(db_3pi.columns)[0]] == idNum).dropna()
            eventRa = event['R.A.'].split(',')[0] #values gives np arrays
            eventDec = event['Dec.'].split(',')[0] #'hh:mm:ss.sss'

        # converting to degrees
        eventCoords = SkyCoord(eventRa, eventDec, unit=(u.hourangle, u.deg))
        eventRa = eventCoords.ra.deg
        eventDec = eventCoords.dec.deg
#TODO: ircs or fk5

        # get event pixel coords
        w = WCS(filename)
        eventX, eventY = w.all_world2pix(eventRa, eventDec, 1)
#TODO check if header for fits is ok
        # to get image dimensions in wcs:
        maxX = image_file[0].header['NAXIS1']
        maxY = image_file[0].header['NAXIS2']
        maxRa, maxDec = w.all_pix2world(1,maxY,1)
        minRa, minDec = w.all_pix2world(maxX,1,1)

        # fix formatting
        maxRa = maxRa.item(0)*u.deg
        minRa = minRa.item(0)*u.deg
        maxDec = maxDec.item(0)*u.deg
        minDec = minDec.item(0)*u.deg
        
        if abs(maxRa/u.deg - eventRa) > 0.01 or \
            abs(minRa/u.deg - eventRa) > 0.01 or \
            abs(maxDec/u.deg - eventDec) > 0.01 or \
            abs(minDec/u.deg - eventDec) > 0.01:
                errorProtocol("Check for coordinate conversion error")



        ''' extract objects '''
        image_data = image_file[0].data

        #fix byte order
        swappedData = image_data.byteswap(True).newbyteorder()
        
        global bkg
        # subtracting out background
        bkg = sep.Background(swappedData)
        if SUBTRACT_BACKGROUND:
            swappedData = swappedData - bkg

        def recursiveExtraction(attemptedThresh):

            global objects
            global bestCandidate
            global magnitude
            global sdssTable
            global chanceCoincidence
            global separation
            global filename
            global blacklist
            global swappedData
            global image_file
            global ra
            global dec
            global hostRa
            global hostDec
            global eventRa
            global eventDec
            global colRealMags
            global colMyMags
            global colFluxes
            global zdict
            global all_myMags
            global all_realMags
            global BAD_IMAGES
            global GOOD_IMAGE
            global m_0s
            global m_0
            global size
            global R_effective
            global hostsData
            global all_myMagsFiltered
            global all_realMagsFiltered
            global all_redshifts
            global all_kronMags
            global all_kronRads
            global kronrad
            global flux
            global z
            global photozs
            global eventz

            global segmap

            objects, segmap = sep.extract(swappedData, attemptedThresh,
                                          err=bkg.globalrms, #var=noise_data2 ,
                                  minarea = MINAREA, deblend_cont = DEBLEND_CONT,
                                  segmentation_map = True)

            # how to calculate kron radius and flux from
            # https://sep.readthedocs.io/en/v1.0.x/apertures.html
            kronrad, _krflag = sep.kron_radius(swappedData, objects['x'], objects['y'],
                                               objects['a'], objects['b'],
                                               objects['theta'], 6.0)

            # objects to be ignored in host detection
            blacklist = set()

            for i in range(len(objects)):
                # remove any 'nan' kronrads and blacklist the objects
                #check if ki is nan:
                if np.isnan(kronrad[i]):
                    #if an object near event (likely host) fails,
                    # redo getting objects using a higher thresh
                    if abs(objects['x'][i] - eventX) < 20 and \
                       abs(objects['y'][i] - eventY) < 20:
                           noteString = "Note: raising threshold to %s %s \n" \
                               %(attemptedThresh + 0.3, filename)
                           print(noteString)
                           with open(ERRORFILE, 'a+') as errorfile:
                               errorfile.write(noteString)
                           recursiveExtraction(attemptedThresh + 0.3)
                           return
                    # otherwise blacklist and give an arbitrary valid kronrad
                    blacklist.add(i)
                    kronrad[i] = 0.1

        ''' '''
        
        # start with default threshold, if it fails then use a higher one
        #three_sigma = 3. * bkg.globalrms

        recursiveExtraction(THRESHOLD)
               
        #remove any nans, replace with cell saturation
        #TODO better solution?
        if not USE_3PI:
            sat = image_file[0].header['HIERARCH CELL.SATURATION']
        else:
#TODO fix
            sat = 100000
        swappedData[np.isnan(swappedData)] = sat
        
        flux = []
        for i in range(len(kronrad)):
            try:
                thisflux, _fluxerr, _flag = sep.sum_ellipse(swappedData, objects['x'][i],
                                        objects['y'][i], objects['a'][i],
                                        objects['b'][i], objects['theta'][i],
                                        2.5*kronrad[i], subpix=1)
                flux.append(thisflux)
            except:

                flux.append(0)
                errorProtocol("Warning: failed flux calculation on " + filename)
        flux = np.array(flux)

        if idNumString not in sdssTableGroupIndex:
            sdssTable = None
        else:
            sdssTable = fullSdssTable.groups[sdssTableGroupIndex[idNumString]]

        # the last argument is 1 b/c first pixel of FITS should be 1,1 not 0,0?
        ra, dec = w.all_pix2world(objects['x'], objects['y'], 1)
        objCoords = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs') # coords of all objs

        # for collecting mags and fluxes to calculate the zero for this file
        colRealMags = []
        #colRealMags = Table([[]]*8, names=magNames)
        colFluxes = []
        photozs = [None]*len(objects)
        for i in range(len(objects)):
            curRa = ra[i]
            curDec = dec[i]
            if sdssTable:
                # for each object, iterate through table until finding object
                #at that location
                for j in range(len(sdssTable)):
                    # check if sdss object matching object i's location
                    if abs(sdssTable['ra'][j] - curRa)*u.deg < MINDIST\
                        and abs(sdssTable['dec'][j] - curDec)*u.deg < MINDIST:

                        # if there's a valid photz, store for future comparison
                        if sdssTable['type'][j] == 3 \
                           and sdssTable['z'][j] \
                           and not np.isnan(sdssTable['z'][j]) \
                           and sdssTable['z'][j] != -9999.0:
                               photozs[i] = sdssTable['z'][j]

                        # if its a star, and its not right on the event, blacklist
                        if sdssTable['type'][j] == 6 \
                            and not abs(sdssTable['ra'][j] - eventRa)*u.deg < MINDIST \
                            and not abs(sdssTable['dec'][j] - eventDec)*u.deg < MINDIST:
                            blacklist.add(i)

                        #if its a star brighter than 21, use its magnitude in zero point calculation
                        if sdssTable['type'][j] == 6 and sdssTable[FILTERS[filterNum]][j] < 21. \
                            and objects['x'][i] - 2.5 * kronrad[i] >= 0 \
                            and objects['x'][i] + 2.5 * kronrad[i] <= maxX \
                            and objects['y'][i] - 2.5 * kronrad[i] >= 0 \
                            and objects['y'][i] + 2.5 * kronrad[i] <= maxY:
                                colRealMags.append(sdssTable[FILTERS[filterNum]][j])
                                colFluxes.append(flux[i])
                        break

        try:
            colFluxes = np.array(colFluxes)
            if not USE_3PI: #convert to flux from counts 
                colFluxes = colFluxes/float(image_file[0].header['EXPTIME'])

#            magcache = [idNum, filterNum]
            m_0s = colRealMags + 2.5 * np.log10(colFluxes)
            m_0s = m_0s[~np.isnan(m_0s)] #remove nans
            if m_0s.any():
                m_0 = np.median(m_0s)
                global m0collector
                m0collector[filterNum].append(m_0)
                
      
            continue
        except:
            raise
    with open("fixed_m0s.txt", 'w+') as f:
        f.write(str(m0collector))
        
    np.save("fixed_m0s", np.array(m0collector), allow_pickle=True)
    return

def main():
#TODO remove
    import time
    start = time.time()
    #figure out which files to use based on value specified at top
    if FILES == 'all' or FILES =='range' or FILES =='new random':
    
        if not USE_3PI:
            filenames = sorted(glob.glob(SOURCEDIR + '/ps1hosts/psc*.[3-6].fits'))
        else:
            filenames = sorted(glob.glob(SOURCEDIR + '/3pi_hosts/ps_*.[3-6].fits'))
            filenames_3pi = sorted(glob.glob(SOURCEDIR + '/3pi_hosts/3pi*.[3-6].fits'))
            
        # If USE_3PI, in addition to usual extraction on orig. data, extract on  3pi additions 
        if FILES == 'range':
            extraction(filenames[RANGE[0]:RANGE[1]])
            if USE_3PI:
                extraction(filenames_3pi[RANGE[0]:RANGE[1]], additional=True)
                
        elif FILES == 'new random':
            randFiles = []
            for i in range(20):
                randFiles.append(filenames[int(len(filenames)*random.random())])
            extraction(randFiles)
            if USE_3PI:
                randFiles_3pi = []
                for i in range(20):
                    randFiles_3pi.append(filenames_3pi[int(len(filenames_3pi)*random.random())])
                extraction(randFiles_3pi, additional=True)
                
        elif FILES =='all':
            extraction(filenames)
            if USE_3PI:
                extraction(filenames_3pi, additional=True)
    elif FILES == 'preset random':
        if USE_3PI: 
            errorProtocol("WARNING PRESET RANDOM OPTION DOES NOT USE 3PI ADDITIONS")
        from presetrandomfiles import fileset
        extraction(fileset[:60])
        
    elif FILES == 'nonsquare':
        from nonsquare import fileset
        extraction(fileset)
    elif FILES == 'specified':
        if SPECIFIED[0][-16:-13] == '3pi':
            addl = True
            if not USE_3PI:
                raise Exception("Specified 3pi image, must set USE_3PI")
        elif SPECIFIED[0][-16:-14] == 'ps':
            addl = False
            if SPECIFIED[0][-16:-13] == 'psc' and USE_3PI:
                raise Exception("Specified NON 3pi image, must unset USE_3PI")
            elif SPECIFIED[0][-16:-13] == 'ps_' and not USE_3PI:
                raise Exception("Specified 3pi image, must set USE_3PI")
        extraction(SPECIFIED, additional=addl)
    else:
        raise Exception('invalid FILE specification')
    end = time.time()
    print(end - start)
    global m0collector
    m0collector = np.array(m0collector)
    np.save("mOcollector", m0collector)

if __name__ == "__main__":
     main()
