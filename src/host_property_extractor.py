# -*- coding: utf-8 -*-

#TODO: all constants at top
#TODO: vertical line
'''
fits files must be named pscxxxxxx.f.fits
and FILTERS must be such that filter number f corresponds with filter filters[f]
alertstable_vs and alertstablevs.lasthalf must be in directory for sn locations
sdssTableGroupIndex.py, sdss_queries.dat

'''

import numpy as np
import pandas as pd
import os
import glob
import random
#import csv
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
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table
#from sdssTableGroupIndex import sdssTableGroupIndex
import json

PROJ_HOME = os.environ['DATA_SRCDIR']

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--mask', action='store_true', 
                    help='only generate mask and date files for cnn')
parser.add_argument('--use_prev', action='store_true', 
                    help='This is not the first time script has been run; \
                        use mag zero points and property averages from last run')
args = parser.parse_args()


if args.mask:
    print("mask making only")
    
if args.use_prev:
    print("Using m_0s and averages from a prev. run")
    
#TODO comment all constants

#TODO fix!!!!
SOURCEDIR = PROJ_HOME + '/src/ps1hosts' #'/mnt/c/Users/Noel/Desktop/summer2019/src/ps1hosts' #pics location
DICTDIR = PROJ_HOME + '/src' #/mnt/c/Users/Noel/Desktop/summer2019/src' #data files location
OUTPUT_DIR = PROJ_HOME + '/src/outputs'
ERRORFILE = OUTPUT_DIR + '/errorfile.txt'

with open(OUTPUT_DIR + '/sdss_queries_index.txt') as f:
    sdssTableGroupIndex = json.load(f)


#os.getcwd() #"/mnt/d/Summer 2019 Astro" 
#"C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd()
FILLER_VAL = None
THRESHOLD = 3 #sigma of detection
MAXTHRESH = 30 # Do not raise threshhold beyond this
PSF = 4 #the FWHM
MINAREA = 3 * (PSF/2)**2
DEBLEND_CONT = 0.01 # for sep.extract. 1.0 to turn off deblending, 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as
FILTERS = [None, None, None, 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z']
LIKELIHOOD_THRESH = 0.2 # only choose hosts with chance coincidence must be below this 
TYPES = ['SNIIn', 'SNIa', 'SNII', 'SNIbc', 'SLSNe']
FILENAME_PREFIX = SOURCEDIR + "/psc"
LOWEST_MAG = 26 #limiting mag is 25
#TODO make note: to change a per-filter property, change in 3 places
#TODO make each property into an object?


#USAGE FLAGS:
#TODO MAKE ARGPARSER!!!!
WRITE_CSV = OUTPUT_DIR + "/galaxiesdata.csv" # filename to write to or None
#MAG_TEST_ALL = False

CHECK_DISTANCE = 5 #print all files with most likely host farther than this arcsecs
PLOT_ALL = False
PLOT_ERR =  True #plots only files that give errors or low probability
PLOT_DIR = OUTPUT_DIR + '/plots' # where to put plot images
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)
ONLY_FLAG_ERRORS = True # catch errors, print filename, move on
FILES = 'all' #options are 'all', 'preset random', 'new random', 'range', 
#'specified', 'nonsquare'

if args.use_prev:
    if os.path.isfile(OUTPUT_DIR + '/m0collector.npy'):
        prev_m0s = np.load(OUTPUT_DIR + '/m0collector.npy', allow_pickle=True)
        FILTER_M0S = [0]*3
        for i in prev_m0s[3:]:
                FILTER_M0S.append(np.median(i))
    else:
        raise Exception("Could not find m0collector.npy from previous run. \
                            Run again without --use_prev flag?")
    if os.path.isfile(WRITE_CSV):
        prev_csv = pd.read_csv(WRITE_CSV)
        AVRG_OFFSETS = [None]*7
        AVRG_RADII = [None]*7
        AVRG_ANGLES = [None]*7
        AVRG_PIXELRANKS = [None]*7
        AVRG_ELLIPTICITIES = [None]*7
        AVRG_REDSHIFTS = [None]*7
        
        good_rows = np.where(prev_csv['host_found_3']==1)
        def get_good_med(header):
            arr = prev_csv[header]
            return np.median(arr[good_rows])
        for i in range(3,7):
            AVRG_OFFSETS[i] = get_good_med('separation (kpc)_%s'%i)
            AVRG_RADII[i] =  get_good_med('KronRad (kpc)_%s'%i)
            AVRG_ANGLES[i] =  get_good_med('Angle_%s'%i)
            AVRG_PIXELRANKS[i] =  get_good_med('pixelRank_%s'%i)
            AVRG_ELLIPTICITIES[i] =  get_good_med('Ellipticity_%s'%i)
            AVRG_REDSHIFTS[i] =  get_good_med('redshift_%s'%i)

    else:
        raise Exception("Count not find %s from previous run. \
                            Run again without --use_prev flag?")
else:
    print("using hardcoded default m0s")
    FILTER_M0s = (None, None, None, 22.918, 22.822, 23.652, 23.540) 
# FILTER_M0s[i] is the average magnitude zero point in filter i, 
#use by default if no stars for calibration
# Also use if star calibration yields outlier (> 1 away from average value)

#TODO delete
SPECIFIED = []
to_check = [160103, 180313, 590123, 50296, 90034, 50601]
for f in to_check:
    SPECIFIED.extend(glob.glob((SOURCEDIR + '/psc*%i*.[3-6].fits' % f)))

SPECIFIED = [SOURCEDIR + '/psc020121.3.fits', 
             SOURCEDIR + '/psc020121.4.fits',
             SOURCEDIR + '/psc020121.5.fits',
             SOURCEDIR + '/psc020121.6.fits']
RANGE = (0,3)
m0collector = [None, None, None, [], [], [], []]
BAD_COUNT = 0
'''make header'''
COLUMNS =['ID', 'hostRa', 'hostDec', 'offby', 'hectoZ', 'redshift_dif']

perImageHeaders = ['KronRad (kpc)', 'separation (kpc)', 'area (kpc^2)', 'sep/sqrt(area) (kpc)',
                   'x', 'y','KronMag', 'Abs. Mag', 'Angle',
                   'Ellipticity', 'RA',  'DEC', 
                   'Discrepency (arcsecs)', 'pixelRank', 'redshift', 'chanceCoincidence',
                   'host_found']
for i in range(3,7):
    for val in perImageHeaders:
        COLUMNS.append(val + '_' + str(i))

# clear any previous logs from errorfile
if os.path.exists(ERRORFILE):
  os.remove(ERRORFILE)
  
  
#all_myMags = []
#all_realMags = []

#for checking how many are in their host
# just out of curiosity
INSIDE = {}
OUTSIDE = {}
INSIDE['SNIa'] = set()
INSIDE['SNIbc'] = set()
INSIDE['SNII'] = set()
INSIDE['SNIIn'] = set()
INSIDE['SLSNe'] = set()

OUTSIDE['SNIa'] = set()
OUTSIDE['SNIbc'] = set()
OUTSIDE['SNII'] = set()
OUTSIDE['SNIIn'] = set()
OUTSIDE['SLSNe'] = set()


# for naming plots files
namecount = 0
def namecountgen():
    global namecount
    namecount += 1
    return namecount

#Four Image objects are created per Supernova object, one for each filter
class Image:    
    def __init__(self, idNumString, filterNum, event):
        self.idNumString = idNumString
        self.idNum = int(idNumString)
        self.filterNum = filterNum
        self.event = event
        self.eventz = float(zdict[self.idNumString])       
        self.bestCandidate = None
        
    # Extract and collect data from image
    def run(self):
        filename = FILENAME_PREFIX + "%s.%s.fits" % (self.idNumString, self.filterNum)
        w = WCS(filename)
        image_file = fits.open(filename)
        
        # to get image corner coordinates in wcs:
        self.maxX = image_file[0].header['NAXIS1']
        self.maxY = image_file[0].header['NAXIS2']
        maxRa, maxDec = w.all_pix2world(1,self.maxY,1)
        minRa, minDec = w.all_pix2world(self.maxX,1,1)
    
#TODO check correctness
        # fix formatting
        maxRa = maxRa.item(0)*u.deg
        minRa = minRa.item(0)*u.deg
        maxDec = maxDec.item(0)*u.deg
        minDec = minDec.item(0)*u.deg
        
        ''' extract objects '''
        image_data = image_file[0].data

        #fix byte order
        self.swappedData = image_data.byteswap(True).newbyteorder()
        # subtracting out background
        self.bkg = sep.Background(self.swappedData)
        if SUBTRACT_BACKGROUND:
            self.swappedData = self.swappedData - self.bkg
            
#TODO does this need to be after the recursive extraction?                
        
        # If there is a bright/oversaturated object too close by to an object,  
        # and extraction threshold is too low, SEP kronrad calculation will fail.
        # Function will then call itself to retry with threshold raised by 0.3.
        # Higher thresh. may unfortunately cause failure to detect fainter objects.
        # as well as lowered detected magnitude and possibly kronradius, 
        # which are currently not corrected for.
        def recursiveExtraction(attemptedThresh):
            # Segmap: Array of integers with same shape as data. Pixels not belonging 
            # to any object have value 0. All pixels belonging to the i-th 
            # object (e.g., objects[i]) have value i+1
            # used for calculating pixel rank 
            # TODO and possibly detecting if event location is in an object
            self.objects, self.segmap = sep.extract(self.swappedData, attemptedThresh,
                                          err=self.bkg.globalrms, 
                                  minarea = MINAREA, deblend_cont = DEBLEND_CONT,
                                  segmentation_map = True)

            # how to calculate kron radius and flux from
            # https://sep.readthedocs.io/en/v1.0.x/apertures.html
            self.kronrad, _krflag = sep.kron_radius(self.swappedData, self.objects['x'], self.objects['y'],
                                               self.objects['a'], self.objects['b'],
                                               self.objects['theta'], 6.0)

            # self.objects to be ignored in host detection
            self.blacklist = set()

            for i in range(len(self.objects)):
#TODO incorrect@!!!!!!!!!!!!!
                # An objcect may have a 'nan' kronrad for 2 reasons. 
                # Usually because object is too long and narrow or cut off, 
                # unlikely to be the host anyway, or an oversaturated object 
                # with nan pixel values. If object is near event, likely host galaxy
                # is getting merged with an oversaturated object, raise 
                # threshhold until they are separated 
                if np.isnan(self.kronrad[i]):
                    if abs(self.objects['x'][i] - self.event['x']) < 20 and \
                       abs(self.objects['y'][i] - self.event['y']) < 20:
                           noteString = "Note: raising threshold to %s \n" \
                               %(attemptedThresh + 0.3)
                           self.errorProtocol(noteString)
                           if attemptedThresh >= MAXTHRESH:
                               return
                           else:
                               recursiveExtraction(attemptedThresh + 0.3)
                               return
                    # otherwise blacklist and give an arbitrary valid kronrad
                    self.blacklist.add(i)
                    self.kronrad[i] = 0.1

        # start with default threshold
        recursiveExtraction(THRESHOLD)
        
        # remove any nan pixels in image, which are due to saturation
        # replace with cell saturation value for closest magnitude accuracy
        # this must be done after recursiveExtraction to make sure we are 
        # separating host from oversaturated objects
        #TODO better solution?
        sat = image_file[0].header['HIERARCH CELL.SATURATION']
        self.swappedData[np.isnan(self.swappedData)] = sat
        
        flux = []
        for i in range(len(self.kronrad)):
            try:
                thisflux, _fluxerr, _flag = self.sum_ellipse_extrapolated(i)
                flux.append(thisflux)
            except:
                # probably object was cut off at edge of image. Ignore object
#TODO make sure ok
#TODO was it flux or kronrad which failed on narrow objects?
                self.blacklist.add(i)
                flux.append(0)
                self.bestCandidate = i
                self.errorProtocol("Warning: failed flux calculation on ")
        flux = np.array(flux)

        # fullSdssTable should contain pre-saved querries of event locations, 
        # with objects grouped by query
        # sdssTableGroupIndex maps a group in the table to an event ID number
        if self.idNumString not in sdssTableGroupIndex:
            # no objects in SDSS near this event location
            sdssTable = None
        else:
            sdssTable = fullSdssTable.groups[sdssTableGroupIndex[self.idNumString]]

#TODO the last argument is 1 b/c first pixel of FITS should be 1,1 not 0,0?
#TODO icrs?
            
        # get WCS coords of all detected objects
        self.ra, self.dec = w.all_pix2world(self.objects['x'], self.objects['y'], 1)
        self.objCoords = SkyCoord(self.ra*u.deg, self.dec*u.deg, frame='icrs') # coords of all objs

        # for collecting mags and fluxes to calculate the zero point for this file
        colRealMags = []
        colFluxes = []
        
#TODO check magnitude correctness
        self.photozs = [None]*len(self.objects)
        self.photozerrs = [None]*len(self.objects)
        for i in range(len(self.objects)):
            curRa = self.ra[i]
            curDec = self.dec[i]
            if sdssTable:
                # for each object, iterate through table until finding object
                #at that location
                for j in range(len(sdssTable)):
                    # check if sdss object matching object i's location
                    if abs(sdssTable['ra'][j] - curRa)*u.deg < MINDIST\
                        and abs(sdssTable['dec'][j] - curDec)*u.deg < MINDIST:

                        # if there's a valid photoz, store for future comparison
                        # TYPE 3 is a galaxy, TYPE 6 is a star
                        if sdssTable['type'][j] == 3 \
                           and sdssTable['z'][j] \
                           and not np.isnan(sdssTable['z'][j]) \
                           and sdssTable['z'][j] != -9999.0:
                               self.photozs[i] = sdssTable['z'][j]
                               self.photozerrs[i] = sdssTable['zErr'][j]

                        # if its a star (coded by type 6), 
                        # and its not right on the event, blacklist
                        if sdssTable['type'][j] == 6 \
                            and not abs(sdssTable['ra'][j] - self.event['ra'])*u.deg < MINDIST \
                            and not abs(sdssTable['dec'][j] - self.event['dec'])*u.deg < MINDIST:
                            self.blacklist.add(i)

                        # if it's a star brighter than 21 
                        # and it's flux range is not cut off by image edges,
                        if sdssTable['type'][j] == 6 and sdssTable[FILTERS[self.filterNum]][j] < 21. \
                            and self.objects['x'][i] - 2.5 * self.kronrad[i] >= 0 \
                            and self.objects['x'][i] + 2.5 * self.kronrad[i] <= self.maxX \
                            and self.objects['y'][i] - 2.5 * self.kronrad[i] >= 0 \
                            and self.objects['y'][i] + 2.5 * self.kronrad[i] <= self.maxY:
                                colRealMags.append(sdssTable[FILTERS[self.filterNum]][j])
                                colFluxes.append(flux[i])
                        break

        colFluxes = np.array(colFluxes)

        # calculate zero points from stars using SDSS magnitudes (colRealMags)
        self.exposure_time = float(image_file[0].header['EXPTIME'])
        m_0s = colRealMags + 2.5 * np.log10(colFluxes/float(self.exposure_time))
        m_0s = m_0s[~np.isnan(m_0s)] # remove nans
        if m_0s.any():
            self.m_0 = np.median(m_0s)
            global m0collector # for calculating median zero points for future runs
            m0collector[self.filterNum].append(self.m_0)

        #no SDSS stars to calculate zero point
        else:
                self.m_0 = FILTER_M0s[self.filterNum]
                colFluxes = np.array([])
                colRealMags = np.array([])

        # if detected m_0 is an outlier, use default instead
        if abs(self.m_0 - FILTER_M0s[self.filterNum]) > 1:
                self.m_0 = FILTER_M0s[self.filterNum]

        # calculate all magnitudes using the new zero point
        self.magnitude = -2.5 * np.log10(flux/self.exposure_time) + self.m_0
        self.magnitude[np.where(np.isnan(self.magnitude))] = LOWEST_MAG

#        #For debugging zero point detection. Can be removed
#        if MAG_TEST_ALL:
#            colMyMags = -2.5 * np.log10(colFluxes/self.exposure_time) + self.m_0
#            all_myMags[self.filterNum].extend(colMyMags[:])
#            all_realMags[self.filterNum].extend(colRealMags[:])
#            for k in range(len(colMyMags)):
#                if abs(colMyMags[k] - colRealMags[k]) > 2:
#                    self.errorProtocol("outlier: ")

        '''Chance Coincidence Calculation'''
        # size is half light radius
        size, _flag = sep.flux_radius(self.swappedData, self.objects['x'], self.objects['y'],
                          6.*self.objects['a'], 0.5, normflux=flux, subpix=5)

        # convert to arcsecs
        degPerPix = image_file[0].header['CDELT1']
        size = size * degPerPix *3600

        self.separation = self.objCoords.separation(self.event['coords']).arcsec
        
        # Observed number density of galaxies brighter than magnitude M (From Berger 2010)
        sigma = 10 ** (0.33 * (self.magnitude - 24) - 2.44) / (0.33 * np.log(10))
        # Effective radius
        R_effective = np.sqrt(np.abs(self.separation)**2 + 4 * np.abs(size) ** 2)

        # Probability of chance coincidence
        self.chanceCoincidence = 1 - np.exp(-np.pi * R_effective ** 2 * sigma)
        # exclude any blacklisted by setting chanceCoincidence to 1, 
        # so will never be chosen as most likely candidate
        for i in self.blacklist:
            self.chanceCoincidence[i] = 1
        
        # Calculate ellipticities
        self.ellipticity = 1 - (self.objects['b']/self.objects['a'])

        # Calculate absolute quantities, convert to kpc, and 
        # remove units for dimensionless math
        dA = cosmo.angular_diameter_distance(self.eventz)*1000/u.Mpc # in kpc.
        self.area = self.objects['npix'] * (degPerPix*(np.pi/180.)*dA)**2 #kpc^2
        self.dL = cosmo.luminosity_distance(self.eventz)*1000/u.Mpc # in kpc
        self.absMag = self.magnitude - 5*np.log10(self.dL) - 10 + 2.5 * np.log10(1.+self.eventz)
        # separation is in arcseconds
        self.separationKpc = self.separation * (1./3600.) * (np.pi/180.)*dA
        self.kronradKpc = self.kronrad * degPerPix*(np.pi/180.)*dA
        self.default_dist = PSF/2*degPerPix*(np.pi/180.)*dA
        self.absLimMag = LOWEST_MAG - 5*np.log10(self.dL) - 10 + 2.5 * np.log10(1.+self.eventz)
        image_file.close()
        
# TODO filter out all filler data before smoting
   
    # returns true IFF pixel corresponding to event location is in object
    def isEventIn(self, objNum):
        x, y = int(self.event['x']), int(self.event['y'])
        # NOTE SEP'S SEGMAP INDEXES Y-FIRST; INDICES ARE ROTATED SUCH THAT 
        # OBJECT N IS CENTERED AROUND segmap[objects['y'][N]][objects['x'][N]]
        return self.segmap[y][x]==objNum + 1

    
    #sets self.bestCandidate. Returns True if chosen by matching redshift,
    # meaning object with lowest chanceCoincidence was not touching event.
    # returns False otherwise
    def attemptBestCandidate(self):
        if len(self.objects) == 0: # no objects were detected in this image
            self.bestCandidate = None
            return False
        
        self.bestCandidate = np.argmin(self.chanceCoincidence)
        photoz_matched = False # whether bestCandidate was chosen by photoz
        # If a candidate not touching event location is chosen, check if any 
        # other objects have photoz matching event redshift, use if so
#TODO ask Ashley if this is okay
        if not self.isEventIn(self.bestCandidate):
        #if self.separation[self.bestCandidate] > CHECK_DISTANCE:
            if self.chanceCoincidence[self.bestCandidate] < 0.02:
                self.errorProtocol("warning, checking for matching photozs for far but likely candidate")
             
            # check for objects with photoz within 10% of event redshift
            # choose matching redshift object with lowest chance coincidence
            for k in range(len(self.photozs)):
                if self.photozs[k] and abs(self.photozs[k] - self.eventz)/self.eventz < 0.1:
                    if photoz_matched:
                        self.errorProtocol("multiple matching photozs")
                        if self.chanceCoincidence[k] < self.chanceCoincidence[self.bestCandidate]:
                               self.bestCandidate = k
                        photoz_matched = True
                  
        return photoz_matched
    
    # Like sep.sum_ellipse but extrapolates if cut-off objects.
    # Assumes symmetry. If cut-off, returns 4*flux in most-intact quarter
    # otherwise, use sep.sum_ellipse
    def sum_ellipse_extrapolated(self, objNum):
        data = self.swappedData
        x0 = self.objects['x'][objNum]
        y0 = self.objects['y'][objNum]
        a = self.objects['a'][objNum]
        b = self.objects['b'][objNum]
        theta_orig = self.objects['theta'][objNum]
        r = 2.5*self.kronrad[objNum]
        
        theta = np.pi/2 - theta_orig
        if theta == 0.:
            theta=0.0001 #to avoid divide by 0 errors
        copy = np.copy(data)
        x_len = len(data)
        y_len = len(data[0])
        
        def in_ellipse(x,y,r):
            return ((x-x0)*np.cos(theta) + (y-y0)*np.sin(theta))**2/a**2 + \
                ((x-x0)*np.sin(theta) - (y-y0)*np.cos(theta))**2/b**2 < r**2 
                
        def in_quarter(x,y,q_key):
            # q_key of 0: lower x side is cutoff, use quarter with higheset x values
            if q_key==0:
                    return ((x-x0) > (y-y0)/np.tan(theta)) and\
                           ((x-x0) > (y-y0)/ -1 * np.tan(theta))
            elif q_key==1: #use lowest x value
                    return ((x-x0) < (y-y0)/np.tan(theta)) and\
                           ((x-x0) < (y-y0)/ -1 * np.tan(theta))
            elif q_key==2: #use highest y vals
                    return (y-y0 > np.tan(theta)*(x-x0)) and\
                           (y-y0 > -1/np.tan(theta)*(x-x0))
            elif q_key==3: # use lowest y vals
                    return (y-y0 < np.tan(theta)*(x-x0)) and\
                           (y-y0 < -1/np.tan(theta)*(x-x0))  
            else:
                raise ValueError("Invalid quartile key")
                                
        #check edges to see where elliipse gets cut off, 
        cutoff_sides = [0]*4 #num. pixels in ellipse crossing [x=0, x=xmax, y=0, y=ymax]
        for y in range(y_len):
            if in_ellipse(0,y,1): #ellipse extends past x=0 edge
                cutoff_sides[0] += 1
            if in_ellipse(x_len-1,y,1):#ellipse extends past x=xmax edge
                cutoff_sides[1]+=1
                
        for x in range(x_len):
            if in_ellipse(x,0,1): #ellipse extends past y=0 edge
                cutoff_sides[2] += 1
            if in_ellipse(x, y_len-1,1):#ellipse extends past y=ymax edge
                cutoff_sides[3]+=1
                
        # Object is not cut off so we can just use sep.sum_ellipse
        if not np.any(cutoff_sides):   
            return sep.sum_ellipse(data,x0, y0, a, b, theta_orig, r, subpix=1)
        else:
            
            # set all pixels outside object flux-ellipse bounds to 0
            for x in range(len(data)):
                for y in range(len(data[0])):
                    if not in_ellipse(x,y,r):
                        copy[x][y] = 0
            
            # set everything outside most-intact-quarter to 0
            cutoff_side = np.argmax(cutoff_sides)
            for x in range(x_len):
                for y in range(y_len):
                    if not in_quarter(x,y,cutoff_side):
                        copy[x][y]=0
                        
            #four times flux in most-intact-quarter of galaxy
            return 4. * np.sum(copy)

    # So all filters will use the same selected host galaxy. goodcoords are the
    # coordinates of the most likely candidate out of all 4 filters. Best 
    # candidate was chosen in used_filter. 
    # Sets bestCandidate to object matching goodcoords, if none return None
    def correct_bestCandidate_to(self, goodcoords, used_filter):
        if used_filter == self.filterNum:
            return
        self.bestCandidate = None
        for k in range(len(self.objCoords)):
            if self.objCoords[k].separation(goodcoords) < MINDIST:
                self.bestCandidate = k
    
    # use data from the filter in which host was chosen, but update magnitude
    # Because no objects in location of chosen host were detected in current filter
    def modifyFinalDataFrom(self, data, newFluxParams):
        newFlux, _fluxerr, _flag = sep.sum_ellipse(*newFluxParams, subpix=1)
        newMagnitude = -2.5 * np.log10(self.exposure_time) + self.m_0
#TODO make sure this is correct: old pixel rank used
        #newPixelRank = self.getPixelRank(target=used_bestCandidate, segmap=used_segmap)
        newAbsMag = newMagnitude - 5*np.log10(self.dL) - 10
        
        f = self.filterNum 
        old_filternum = list(data.keys())[0][-1]
        oldMagnitude = data['KronMag_%s' % old_filternum]
        new_data = {}
        for (k, v) in data.items():
            new_key = k[:-1] + str(self.filterNum)
            new_data[new_key] = v
        new_data['KronMag_%s' % f] = newMagnitude
        new_data['Abs. Mag_%s' % f] = newAbsMag
        #new_data['pixelRank_%s' % f] = newPixelRank
        new_data['KronRad (kpc)_%s' % f] = newMagnitude
        oldChanceCoincidence = data['chanceCoincidence_%s' % old_filternum]
#TODO check length, make sure this works
        #adjust chanceCoincidence for the new magnitude
        new_data['chanceCoincidence_%s' % f] = 1 - (1 - oldChanceCoincidence)**(10**(0.33*(newMagnitude - oldMagnitude)))
#TODO make sure above math is correct
        new_data['host_found_%s' % f] = 0
        return new_data
    
    # returns host data from this filter to be entered into csv
    def getFinalData(self):    
        finalDict = {}
        f = self.filterNum               
        bestCandidate = self.bestCandidate
        finalDict['KronRad (kpc)_%s' % f] = self.kronradKpc[bestCandidate]
        finalDict['separation (kpc)_%s' % f] = self.separationKpc[bestCandidate]
        finalDict['area (kpc^2)_%s' % f] = self.area[bestCandidate]
        finalDict['sep/sqrt(area) (kpc)_%s' % f] = self.separationKpc[bestCandidate] /\
            np.sqrt(self.area[bestCandidate])
        finalDict['x_%s' % f] = self.objects['x'][bestCandidate] - self.event['x']
        finalDict['y_%s' % f] = self.objects['y'][bestCandidate] - self.event['y']
        finalDict['KronMag_%s' % f] = self.magnitude[bestCandidate]
        finalDict['Abs. Mag_%s' % f] = self.absMag[bestCandidate]
        finalDict['Angle_%s' % f] = self.objects['theta'][bestCandidate]
        finalDict['Ellipticity_%s' % f] =  self.ellipticity[bestCandidate]
        finalDict['RA_%s' % f] = self.ra[bestCandidate]
        finalDict['DEC_%s' % f] = self.dec[bestCandidate]
        finalDict['Discrepency (arcsecs)_%s' % f] = self.dec[bestCandidate]
        finalDict['pixelRank_%s' % f] = self.getPixelRank()
        finalDict['chanceCoincidence_%s' % f] = self.chanceCoincidence[bestCandidate]
        finalDict['host_found_%s' % f] = 1
        finalDict['redshift_%s' % f] = self.eventz
        return finalDict

    # If no likely hosts were detected in any filter for this event
    def getDefaultData(self):
        defaultFinalProperties = {}
        for f in range(3,7):
            defaultFinalProperties['KronRad (kpc)_%s' % f] = self.default_dist #PSF, 2 pixels
            if args.use_prev:
                defaultFinalProperties['separation (kpc)_%s' % f] = \
                    self.default_dist*AVRG_OFFSETS[self.filterNum]/AVRG_RADII[self.filterNum] 
            else: # ideally data from a use_prev run will be used in the end
                # so the following won't really matter
                defaultFinalProperties['separation (kpc)_%s' % f] = self.default_dist/2 # 1 pixel
            defaultFinalProperties['area (kpc^2)_%s' % f] = \
                np.pi * (self.default_dist) ** 2 
            defaultFinalProperties['sep/sqrt(area) (kpc)_%s' % f] = \
                defaultFinalProperties['separation (kpc)_%s' % f] /\
                np.sqrt(defaultFinalProperties['area (kpc^2)_%s' % f])
            defaultFinalProperties['x_%s' % f] = 0
            defaultFinalProperties['y_%s' % f] = 0
            defaultFinalProperties['KronMag_%s' % f] = LOWEST_MAG
            
#TODO             FIXXX to 2sigma
            
            defaultFinalProperties['Abs. Mag_%s' % f] = self.absLimMag
            if args.use_prev:
                defaultFinalProperties['Angle_%s' % f] = AVRG_ANGLES[self.filterNum]
                defaultFinalProperties['Ellipticity_%s' % f] =  AVRG_ELLIPTICITIES[self.filterNum]
            else:
                defaultFinalProperties['Angle_%s' % f] = 0
                defaultFinalProperties['Ellipticity_%s' % f] =  0.3
            defaultFinalProperties['RA_%s' % f] = self.event['ra']
            defaultFinalProperties['DEC_%s' % f] = self.event['dec']
            defaultFinalProperties['Discrepency (arcsecs)_%s' % f] = None
            if args.use_prev:                
                defaultFinalProperties['pixelRank_%s' % f] = AVRG_PIXELRANKS[self.filterNum]
            else:
                defaultFinalProperties['pixelRank_%s' % f] = 0.5
            defaultFinalProperties['chanceCoincidence_%s' % f] = 1
            defaultFinalProperties['host_found_%s' % f] = 0
            if args.use_prev:
                defaultFinalProperties['redshift_%s' % f] = AVRG_REDSHIFTS[self.filterNum]
            else:
                defaultFinalProperties['redshift_%s' % f] = 0.2
        return defaultFinalProperties
    
    # returns |Redshift_SN - Photoz_galaxy| / Photoz error on galaxy 
    # or -1 if no SDSS photoz best candidate galaxy
    def getPhotozDif(self):
        photoz = self.photozs[self.bestCandidate]
        photozerr = self.photozerrs[self.bestCandidate]
        eventz = self.eventz
        
        #PHOTOZ failed or no match
        if photoz == None or photozerr==None or photozerr < -99  or np.isnan(photozerr): 
            return -1
        else:
            return abs(eventz - photoz)/photozerr
        

    # saves and displays image with SDSS identified stars and 
    # blacklisted (disqualified) objects circled in blue,
    # best host candidate circled in green, all other detected objects circled
    # in red, event location marked with purple triangle. Saves to PLOT_DIR
    def plot(self, myVmin=None, myVmax=None, target=None, title=''):
        green = [target] if target else [self.bestCandidate]
        # make the destination directory if it does not exist
        if not os.path.isdir(PLOT_DIR):
            os.mkdir(PLOT_DIR)
        to_circle = range(len(self.objects))
    
        fig, ax = plt.subplots()
    
        if myVmax == None or myVmin == None:
            _im = ax.imshow(self.swappedData, interpolation='nearest', cmap='gray_r')
        else:
            _im = ax.imshow(self.swappedData, interpolation='nearest', cmap='gray_r',
                            vmin = myVmin, vmax = myVmax)
        # triangle on event location
        p = RegularPolygon((int(self.event['x']), int(self.event['y'])), 3, radius=8)
        p.set_edgecolor('purple')
        p.set_facecolor('purple')
        ax.add_artist(p)
    
        # plot an ellipse for each object
        for i in to_circle:
            e = Ellipse(xy=(self.objects['x'][i], self.objects['y'][i]),
                        width=6*self.objects['a'][i],
                        height=6*self.objects['b'][i],
                        lw = 3,
                        angle=self.objects['theta'][i] * 180. / np.pi)
    # TODO check if above angle conversion is correct. Where from?
            e.set_facecolor('none')
            if i in green:
                e.set_edgecolor('green')
            elif i in self.blacklist:
                e.set_edgecolor('blue')
            else:
                e.set_edgecolor('red')
            ax.add_artist(e)
        plt.title(title)
        plt.savefig(PLOT_DIR + "/galaxyimage" + self.idNumString + '_' \
                    + str(namecountgen()) + ".png", dpi=150)
        plt.show()
        plt.close()

    def getPixelRank(self, target=None, segmap=None):
        # NOTE SEP'S SEGMAP INDEXES Y-FIRST; INDICES ARE ROTATED SUCH THAT 
        # OBJECT N IS CENTERED AROUND segmap[objects['y'][N]][objects['x'][N]]
        if not target:
            target = self.bestCandidate
        if not segmap:
            segmap = self.segmap
        a = np.where(segmap == target+1)
        pixels = []
        for k in range(len(a[0])):
            pixels.append((a[0][k], a[1][k])) #list of tuples of coords
            global INSIDE
            INSIDE[typeDict[self.idNumString]].add(self.idNumString)
        if not (int(self.event['y']), int(self.event['x'])) in pixels:
            global OUTSIDE
            OUTSIDE[typeDict[self.idNumString]].add(self.idNumString)
            return 0
        
        def sortkey(x):
            a, b = x
            return self.swappedData[a][b]        
        pixels.sort(key = sortkey)
        location = pixels.index((int(self.event['y']), int(self.event['x'])))
        return float(location)/float(len(pixels))
    
    # loggging
    def errorProtocol(self, e, target=None):
            curFile = self.idNumString + '.' + str(self.filterNum)
            errorString = str(namecount) + '\n' + str(e) + " " + curFile + '\n'
            print(errorString)
            with open(ERRORFILE, 'a+') as errorfile:
                errorfile.write(errorString)

            chosen = target if target else self.bestCandidate
            e = str(e)
            if PLOT_ERR:
                padded_e = e + '                                                                         '
                my_title = padded_e[:30] + '\n' + curFile[-16:]
                self.plot(myVmin = 0, myVmax = 3000, target=chosen, title=my_title)
                self.plot()
            if ONLY_FLAG_ERRORS or e=='far' or e=='unlikely':
                return
            else:
                print("raising")
                raise

def pre_load_dictionaries():
    '''load event type dictionary'''
    global typeDict
    typeDict = {}
    typefile = open(DICTDIR + '/ps1confirmed_added.txt', 'r')
    typefile.readline() #get rid of header
    for line in typefile:
        parts = line.split()
        eventId = parts[0][3:]
        eventType = parts[1]
        typeDict[eventId] = eventType
    typefile.close()
    
    
    ''' load event redshifts' dictionary '''
    global zdict
    zdict = {}
    zfile = open(DICTDIR + '/new_ps1z.dat', 'r')
    zfile.readline() #get rid of heade
    
    for line in zfile:
        parts = line.split()
        eventId = parts[0][3:]
        redshift = float(parts[1])
        zdict[eventId] = redshift
    
    zfile.close()
    #not in zdict, looked up from ps1 website 
    #http://telescopes.rc.fas.harvard.edu/ps1/2014/alerts/
#     zdict['100014'] = 0.357
#     zdict['300220'] = 0.094
#     zdict['380108'] = 0.159
#     zdict['080079'] = 0.151
#     zdict['460103'] = 0.32
#TODO Ask ashley if 460103 is PS1-12c
    
    
    
    '''load event location dictionary'''
    global db
    db = pd.read_csv(DICTDIR + '/alertstable_v3',sep=None,index_col = False,
                   engine='python')
    db2 = pd.read_csv(DICTDIR + '/alertstable_v3.lasthalf',sep=None,index_col = False,
                    engine='python')
    db = db.append(db2,ignore_index=True)
    
    
    
    '''load prerun sdss queries'''
    global fullSdssTable
    fullSdssTable = Table.read(OUTPUT_DIR + '/sdss_queries.dat', format='ascii')
    fullSdssTable = fullSdssTable.group_by('idnum')
    
    
    '''load real host data'''
    global hostsData
    with open(DICTDIR + '/hosts.dat', 'r') as hostsfile:
            hostsData = hostsfile.read()
# not sure if I should be able to json.load this, but it fails and this works
    hostsData = ast.literal_eval(hostsData)

class Supernova:
    
    def __init__(self, idNumString):
        self.idNumString = idNumString
        self.idNum = int(idNumString)
        
    def getSnFinalData(self, chosen_loc):
        finalDict = {'ID':self.idNum}
        hostRa = hostsData[self.idNumString]['host_ra']
        hostDec = hostsData[self.idNumString]['host_dec']
        try:
            hostCoords = SkyCoord(hostRa, hostDec, unit=(u.hourangle, u.deg))
            #convert to decimal deg:
            finalDict['hostRa'] = hostCoords.ra.deg
            finalDict['hostDec'] = hostCoords.dec.deg
            if hostRa and chosen_loc:
                finalDict['offby'] = hostCoords.separation(chosen_loc).arcsec
            else:
                finalDict['offby'] =None
            finalDict['hectoZ'] = hostsData[self.idNumString]['redshift']
        except ValueError:
            if len(hostRa) > 12: # hecto gave multiple host galaxies
                finalDict['hostRa'] = "multiple"
                finalDict['hostDec'] = "multiple"
                finalDict['offby'] = None
                finalDict['hectoZ'] = "multiple"
            else:
                raise
                
#TODO note that redshift_offset of -1 indicates no photoz, positive is difference
        if self.used_default:
            finalDict['redshift_dif'] = -1
        else:
            finalDict['redshift_dif'] = self.images[self.filter_to_use].getPhotozDif()
  
        return finalDict
        
    def run(self):
        # find pixel coordinates of event
        # if image is cutoff and so not square, event won't be in center
        e = db.where(db['eventID'] == self.idNum).dropna()
        eRa = e['ra'].values[0] #values gives np arrays
        eDec = e['dec'].values[0] #'hh:mm:ss.sss'
        #TODO: ircs or fk5
        event = {}
        event['coords'] = SkyCoord(eRa, eDec, unit=(u.hourangle, u.deg))
        event['ra'] = event['coords'].ra.deg
        event['dec'] = event['coords'].dec.deg
        # get event pixel coords
        filename = FILENAME_PREFIX + "%s.3.fits" % self.idNumString
        w = WCS(filename)
        event['x'], event['y'] = w.all_world2pix(event['ra'], event['dec'], 1)
        self.event = event
              
        self.images = [None]*7
        good_images = []
        BAD_IMAGES = []
        good_photozs = []
        for x in range(3,7):
            self.images[x] = Image(self.idNumString, x, event)
            self.images[x].run()
            photozMatched = self.images[x].attemptBestCandidate()
#TODO switch to absolute distance? take into account size? make sure we aren't correcting big galaxies
#TODO mindist vs. checkdist
#TODO What if multiple galaxies are within mindist and dif ones are selected in dif filters?
# I thought I fixed this to check if event was inside!
            if self.images[x].bestCandidate != None and \
                self.images[x].objCoords[self.images[x].bestCandidate].separation(event['coords']) < MINDIST:
                good_images.append(x)
            elif photozMatched:
                good_photozs.append(x)
            else:
                BAD_IMAGES.append(x)

        
        if args.mask:
            global masks
        
        '''choose candidates'''        
        if good_images:
            self.filter_to_use = good_images[0]
#TODO remove
            goodCoords = self.images[self.filter_to_use].\
                objCoords[self.images[self.filter_to_use].bestCandidate]
#TODO do I use object at event coords or good coords?
            for x in good_photozs + BAD_IMAGES:
                self.images[x].correct_bestCandidate_to(goodCoords, self.filter_to_use)
        elif good_photozs:
#TODO pick min chance coincidence !!!
            self.filter_to_use = good_photozs[0]
            goodCoords = self.images[self.filter_to_use].\
                objCoords[self.images[self.filter_to_use].bestCandidate]
            #use object at that location in all images
            for x in range(3,7):
                self.images[x].correct_bestCandidate_to(goodCoords, self.filter_to_use)
#TODO test all combinations of good images, bad_images, no objects images            
        else:
            # find minimum chance coincidence object
            first = self.images[BAD_IMAGES[0]]
            if len(first.objects) > 0:
                minChanceCoincidence = first.chanceCoincidence[first.bestCandidate]
            else:
                # no objects were detected in image
                minChanceCoincidence = 2
            minOwner = BAD_IMAGES[0]
            for num in BAD_IMAGES:
                im = self.images[num]
                if len(im.objects) > 0:
                    bestChance = im.chanceCoincidence[im.bestCandidate]
                else: # no objects were detected in image
                    bestChance = 2
                if bestChance < minChanceCoincidence:
                   minChanceCoincidence = bestChance
                   minOwner = num
            self.filter_to_use=minOwner
            if minChanceCoincidence <= 0.02:
                self.images[self.filter_to_use].errorProtocol(
                        "USING BEST CHANCE: %s" % minChanceCoincidence)
                goodCoords = self.images[self.filter_to_use].\
                    objCoords[self.images[self.filter_to_use].bestCandidate]
                for i in range(3,7):
                    self.images[i].correct_bestCandidate_to(goodCoords, self.filter_to_use)
            else:
                #no good candidates at all, GET DEFAULT DATA AND RETURN
                self.images[self.filter_to_use].errorProtocol(
                        "Best.%s.NO_GOOD_CANIDIDATE" % minChanceCoincidence)
                global BAD_COUNT
                BAD_COUNT += 1
                all_sn_data = {}
                for x in range(3,7):
                    all_sn_data.update(self.images[x].getDefaultData())
#TODO make sure used_default is correct                    
                self.used_default=False
                #get non filter dependent data
                all_sn_data.update(self.getSnFinalData(None))
                print('ID: %s' % self.idNum)
                if args.mask:
                    #global masks
                    masks[self.idNumString] = np.zeros((240,240))
                    #np.save('masks/mask_%s'%self.idNumString, np.zeros((240,240)))
                    print('bad mask %s' % self.idNumString)
                    return#((np.zeros((240,240)), 
                            #[self.images[0], self.images[1], self.images[2], self.images[3]]))
                else:
                    return all_sn_data
            
            
            
        self.used_default = True
        all_sn_data = {}
        chosen_im = self.images[self.filter_to_use]
        # get final data from the known good filter first
        best_data = chosen_im.getFinalData()
        all_sn_data.update(best_data)
        for i in range(3,7):
            if i == self.filter_to_use: 
                if args.mask: # only collecting mask
                    mask = self.images[i].segmap
                    mask = np.where(mask == self.images[i].bestCandidate + 1, 1 - self.images[i].chanceCoincidence[self.images[i].bestCandidate], 0)
                    #global masks
                    masks[self.idNumString] = mask
                    print('good mask %s' % self.idNumString)
                    #np.save('masks/mask_%s'%self.idNumString, mask)
                    return #(mask, 
                            #[self.images[0], self.images[1], self.images[2], self.images[3]])
                else:
                    #already collected data for this filter above
                    pass
            elif self.images[i].bestCandidate == None:
                # the chosen host candidate could not be found in this image
                # modify and use data from filter_to_use which chose the host
#TODO check correctnesss
                usedbestCandidate = self.images[self.filter_to_use].bestCandidate
                newFluxParams =  (self.images[i].swappedData,
                   chosen_im.objects['x'][usedbestCandidate],
                   chosen_im.objects['y'][usedbestCandidate],
                   chosen_im.objects['a'][usedbestCandidate],
                   chosen_im.objects['b'][usedbestCandidate],
                   chosen_im.objects['theta'][usedbestCandidate], 
                   2.5*chosen_im.kronradKpc[usedbestCandidate])
                all_sn_data.update(
                        self.images[i].modifyFinalDataFrom(best_data, newFluxParams))
            else:
                #get filter dependent data
                all_sn_data.update(self.images[i].getFinalData())
                
#TODO make flattening better                
        #get non-filter-dependent data
        chosen_loc = chosen_im.objCoords[chosen_im.bestCandidate]
        all_sn_data.update(self.getSnFinalData(chosen_loc))
        return all_sn_data
    
        if PLOT_ALL:
            for image in self.images[3:]:
                image.plot()
                image.plot(myVmin = 0, myVmax = 1000)        


def extraction(filenames):
    # Sn data tables neccessary for running
    pre_load_dictionaries()
    
    all_redshifts = {}
    all_kronMags = {}
    all_kronRads = {}
    for snType in TYPES:
        all_redshifts[snType] = []
        all_kronMags[snType] = []
        all_kronRads[snType] = []

    all_all_data = []
    
    if args.mask:
        global masks
        masks = {}
        #if not os.path.isdir('masks'):
        #    os.mkdir('masks')
        #masks = [[],[],[],[],[]]
        #dats = [[],[],[],[],[]]
        #numDict = {'SNIa':0, 'SNIbc':1, 'SNII': 2, 'SNIIn':3, 'SLSNe':4}
    for filename in filenames: 
        # filename will be the #3 (g) filter file of the sn
        # check that all 4 filters are present, otherwise skip
        # filename will end with xxx123456.n.fits where n is filternum
        if len(glob.glob(filename[:-6] + '[3-6]' + filename[-5:])) != 4:
            continue
        # to extract transient and filter number from filename of the form
        # */psc190369.6.fits
        dotSplit = filename.split('.')
        # transient identifier is everything after the last 'c' but before
        # the second to last dot
        idNumString = dotSplit[-3].split('c')[-1] #used for getting hectospec data
        # filter number is between second to last dot and last dot
        s = Supernova(idNumString)
#TODO fix idNumString getting
        
        if args.mask:
            print(filename)
            s.run()
#            mask, dat = s.run()
#            dat = np.array(dat)
#            mask = np.array(mask)
#            typ = numDict[typeDict[idNumString]]
#            masks[typ].append(mask)
#            dats[typ].append(dat)
        else:            
            all_all_data.append(s.run())
            
            #TODO column order
            #print(all_all_data)
            df = pd.DataFrame(all_all_data, columns=COLUMNS)
            if WRITE_CSV:
                df.to_csv(WRITE_CSV)
                
    if args.mask:
        np.savez('all_masks', **masks)
        print('saved all masks')
#        for i in range(5):
#            np.save('masks_%s'%i, masks[i])
#            np.save('x_all2_%s'%i, dats[i])
            

def main():
    import time
    start = time.time()
    #figure out which files to use based on value specified at top
    if FILES == 'all' or FILES =='range' or FILES =='new random':
        #CHANGED now only looks for .3 files, assuming the rest are there
        filenames = sorted(glob.glob(SOURCEDIR + '/psc*.3.fits'))
        if FILES == 'range':
            extraction(filenames[RANGE[0]:RANGE[1]])
        elif FILES == 'new random':
            randFiles = []
            for i in range(20):
                randFiles.append(filenames[int(len(filenames)*random.random())])
            extraction(randFiles)
        elif FILES =='all':
            extraction(filenames)
    elif FILES == 'preset random':
        from presetrandomfiles import fileset
        extraction(fileset[:60])
    elif FILES == 'nonsquare':
        from nonsquare import fileset
        extraction(fileset)
    elif FILES == 'specified':
        extraction(SPECIFIED)
    else:
        raise Exception('invalid FILE specification')
    end = time.time()
    print("Time: %s " % (end - start))
    print("BAD COUNT: %s" % BAD_COUNT)
    np.save("m0collector", m0collector)
if __name__ == "__main__":
     main()
     
