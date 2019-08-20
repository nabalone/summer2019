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
from sdssTableGroupIndex import sdssTableGroupIndex

#TODO comment all constants
ERRORFILE = 'errorfile.txt'
SOURCEDIR = os.getcwd() #"/mnt/d/Summer 2019 Astro" #"C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd()
FILLER_VAL = None
THRESHOLD = 3
PSF = 4 #the FWHM
MINAREA = 3 * (PSF/2)**2
DEBLEND_CONT = 0.01 # for sep.extract. 1.0 to turn off deblending, 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as
FILTERS = [None, None, None, 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z']
LIKELIHOOD_THRESH = 0.2
TYPES = ['SNIIn', 'SNIa', 'SNII', 'SNIbc', 'SLSNe']
TYPE_COLORS = {'SNIIn':'co', 'SNIa':'ro', 'SNII': 'bo', 'SNIbc':'go', 'SLSNe': 'mo'}
FILENAME_PREFIX = SOURCEDIR + "/ps1hosts/psc"
LOWEST_MAG = 26 #limiting mag is 25
#TODO make note: to change a per-filter property, change in 3 places
#TODO make each property into an object?

FILTER_M0s = (None, None, None, 22.918, 22.822, 23.652, 23.540) 
# FILTER_M0s[i] is the average magnitude zero point in filter i, use by default if no stars for calibration
# Also use if star calibration yields outlier (> 1 away from average value)

#USAGE FLAGS:
WRITE_CSV = "/galaxiesdata.csv" # filename to write to or None
MAG_TEST_ALL = False
MAG_TEST_STDEV = False
PLOT_REDSHIFTS = False

PRINT_DATA = False

CHECK_DISTANCE = 5 #print all files with most likely host farther than this arcsecs
PLOT_ALL = False
PLOT_ERR =  True #plots only files that give errors or low probability
PLOT_DIR = os.getcwd() + '/plots' # where to put plot images
ONLY_FLAG_ERRORS = True # catch errors, print filename, move on
FILES = 'range' #options are 'all', 'preset random', 'new random', 'range', 'specified', 'nonsquare

#TODO delete
SPECIFIED = []
to_check = [160103, 180313, 590123, 50296, 90034, 50601]
for f in to_check:
    SPECIFIED.extend(glob.glob((SOURCEDIR + '/ps1hosts/psc*%i*.[3-6].fits' % f)))

SPECIFIED = [SOURCEDIR + '/ps1hosts/psc480552.6.fits']
RANGE = (0,100)
m0collector = [None, None, None, [], [], [], []]
BAD_COUNT = 0
'''make header'''
COLUMNS =['ID', 'hostRa', 'hostDec', 'offby', 'hectoZ', 'redshift_dif']

perImageHeaders = ['KronRad (kpc)', 'separation (kpc)', 'area (kpc^2)', 'sep/sqrt(area) (kpc)',
                   'x', 'y','KronMag', 'Abs. Mag', 'Angle',
                   'Ellipticity', 'RA',  'DEC', 
                   'Discrepency (arcsecs)', 'pixelRank', 'chanceCoincidence',
                   'host_found']
for i in range(3,7):
    for val in perImageHeaders:
        COLUMNS.append(val + '_' + str(i))

# remove error file
if os.path.exists(ERRORFILE):
  os.remove(ERRORFILE)
  
  
all_myMags = []
all_realMags = []

# for naming plots files
namecount = 0
def namegen():
    global namecount
    namecount += 1
    return PLOT_DIR + "\galaxyimage" + str(namecount) + ".png"

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
        # remove any nan pixels in image, which are due to saturation
        # replace with cell saturation value for closest magnitude accuracy
        #TODO better solution?
        sat = image_file[0].header['HIERARCH CELL.SATURATION']
        self.swappedData[np.isnan(self.swappedData)] = sat
        
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
                # Usually because object is too long and narrow, 
                # unlikely to be the host anyway.
                # Ocassionally because a nearby bright star disrupted calculation
                # We assume the former and blacklist (ignore) object unless it
                # is near event site, in which case we raise threshold just in case 
                # of the latter
                if np.isnan(self.kronrad[i]):
                    if abs(self.objects['x'][i] - self.event['x']) < 20 and \
                       abs(self.objects['y'][i] - self.event['y']) < 20:
                           noteString = "Note: raising threshold to %s \n" \
                               %(attemptedThresh + 0.3)
                           self.errorProtocol(noteString)
                           recursiveExtraction(attemptedThresh + 0.3)
                           return
                    # otherwise blacklist and give an arbitrary valid kronrad
                    self.blacklist.add(i)
                    self.kronrad[i] = 0.1

        # start with default threshold, if it fails then use a higher one
        recursiveExtraction(THRESHOLD)
        
        flux = []
        for i in range(len(self.kronrad)):
            try:
                thisflux, _fluxerr, _flag = sep.sum_ellipse(self.swappedData, self.objects['x'][i],
                                        self.objects['y'][i], self.objects['a'][i],
                                        self.objects['b'][i], self.objects['theta'][i],
                                        2.5*self.kronrad[i], subpix=1)
                flux.append(thisflux)
            except:
                # probably because object was too narrow. Ignore object
#TODO make sure ok
                self.blacklist.add(i)
                flux.append(0)
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

        #For debugging zero point detection. Can be removed
        if MAG_TEST_ALL:
            colMyMags = -2.5 * np.log10(colFluxes/self.exposure_time) + self.m_0
            all_myMags[self.filterNum].extend(colMyMags[:])
            all_realMags[self.filterNum].extend(colRealMags[:])
            for k in range(len(colMyMags)):
                if abs(colMyMags[k] - colRealMags[k]) > 2:
                    self.errorProtocol("outlier: ")

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
        objMap = np.where(self.segmap == self.objects[objNum]+1)
        return (int(self.event['x']), int(self.event['y'])) in objMap
    
    #sets self.bestCandidate. Returns True if chosen by matching redshift,
    # meaning object with lowest chanceCoincidence was not touching event.
    # returns False otherwise
    def attemptBestCandidate(self):
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

        #check if best candidate is cut off by image boundaries
        if self.objects['x'][self.bestCandidate] + 2.5* self.kronradKpc[self.bestCandidate] > self.maxX or \
            self.objects['x'][self.bestCandidate] - 2.5* self.kronradKpc[self.bestCandidate] < 0 or \
            self.objects['y'][self.bestCandidate] + 2.5* self.kronradKpc[self.bestCandidate] > self.maxY or \
            self.objects['y'][self.bestCandidate] - 2.5* self.kronradKpc[self.bestCandidate] < 0:
                self.errorProtocol("cut off object")                        
#TODO extrapolate                     
        return photoz_matched
    
    def correct_bestCandidate_to(self, goodcoords, used_filter):
        if used_filter == self.filterNum:
            return
        self.bestCandidate = None
        for k in range(len(self.objCoords)):
            if self.objCoords[k].separation(goodcoords) < MINDIST:
                self.bestCandidate = k
#TODO return true or false, check
    
    def modifyFinalDataFrom(self, data, newFluxParams):#, used_bestCandidate, used_segmap):
        # no objects in correct location detected in bad image, 
        #use data from good image but update magnitude
        newFlux, _fluxerr, _flag = sep.sum_ellipse(*newFluxParams, subpix=1)
        newMagnitude = -2.5 * np.log10(self.exposure_time) + self.m_0
#TODO make sure this is correct just save the old pixel rank
        #newPixelRank = self.getPixelRank(target=used_bestCandidate, segmap=used_segmap)
        newAbsMag = newMagnitude - 5*np.log10(self.dL) - 10
        
        f = self.filterNum 
        old_filternum = data.keys()[0][-1]
        oldMagnitude = data['KronMag_%s' % old_filternum]
        new_data = {}
        for (k, v) in data.iteritems():
            new_key = k[:-1] + str(self.filterNum)
            new_data[new_key] = v
        new_data['KronMag_%s' % f] = newMagnitude
        new_data['Abs. Mag_%s' % f] = newAbsMag
        #new_data['pixelRank_%s' % f] = newPixelRank
        new_data['KronRad (kpc)_%s' % f] = newMagnitude
        oldChanceCoincidence = data['chanceCoincidence_%s' % old_filternum]
#TODO check length, make sure this works
        #adjust chanceCoincidence for change in magnitude
        new_data['chanceCoincidence_%s' % f] = 1 - (1 - oldChanceCoincidence)**(10**(0.33*(newMagnitude - oldMagnitude)))
#TODO make sure above math is correct
        new_data['host_found_%s' % f] = 0
        return new_data
    
    def getFinalData(self):    
        '''Get "real" host location and redshifts'''
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
        return finalDict

    def getDefaultData(self):
        defaultFinalProperties = {}
        for f in range(3,7):
            defaultFinalProperties['KronRad (kpc)_%s' % f] = self.default_dist #PSF/2 pixels
            defaultFinalProperties['separation (kpc)_%s' % f] = self.default_dist/2 #1 pixel
            defaultFinalProperties['area (kpc^2)_%s' % f] = \
                np.pi * (self.default_dist) ** 2 #psf/2
            defaultFinalProperties['sep/sqrt(area) (kpc)_%s' % f] = \
                defaultFinalProperties['separation (kpc)_%s' % f] /\
                np.sqrt(defaultFinalProperties['area (kpc^2)_%s' % f])
            defaultFinalProperties['x_%s' % f] = 0
            defaultFinalProperties['y_%s' % f] = 0
            defaultFinalProperties['KronMag_%s' % f] = LOWEST_MAG
            defaultFinalProperties['Abs. Mag_%s' % f] = self.absLimMag
            defaultFinalProperties['Angle_%s' % f] = 0
            defaultFinalProperties['Ellipticity_%s' % f] =  0.3
            defaultFinalProperties['RA_%s' % f] = self.event['ra']
            defaultFinalProperties['DEC_%s' % f] = self.event['dec']
            defaultFinalProperties['Discrepency (arcsecs)_%s' % f] = None
            defaultFinalProperties['pixelRank_%s' % f] = 0.5
            defaultFinalProperties['chanceCoincidence_%s' % f] = 1
            defaultFinalProperties['host_found_%s' % f] = 0
        return defaultFinalProperties
    
    def getPhotozDif(self):
        photoz = self.photozs[self.bestCandidate]
        photozerr = self.photozerrs[self.bestCandidate]
        eventz = self.eventz
        
        #TODO note that redshiftdif of -1 indicates no photoz, positive is difference
        #PHOTOZ failed or no match
        if photoz == None or photozerr==None or photozerr < -99  or np.isnan(photozerr): 
            return -1
        else:
            return abs(eventz - photoz)/photozerr
        

   
    def plot(self, myVmin=None, myVmax=None, target=None):
        green = [target] if target else [self.bestCandidate]
        # make the destination directory if it does not exist
        if not os.path.isdir(PLOT_DIR):
            os.mkdir(PLOT_DIR)
        to_circle = range(len(self.objects))
    
        fig, ax = plt.subplots()
    
        if myVmax == None or myVmin == None:
            _im = ax.imshow(self.swappedData, interpolation='nearest', cmap='gray')
        else:
            _im = ax.imshow(self.swappedData, interpolation='nearest', cmap='gray',
                            vmin = myVmin, vmax = myVmax)
        # triangle on event location
        p = RegularPolygon((int(self.event['x']), int(self.event['y'])), 3, radius=3)
        p.set_edgecolor('purple')
        p.set_facecolor('purple')
        ax.add_artist(p)
    
        # plot an ellipse for each object
        for i in to_circle:
            e = Ellipse(xy=(self.objects['x'][i], self.objects['y'][i]),
                        width=6*self.objects['a'][i],
                        height=6*self.objects['b'][i],
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
        plt.savefig(namegen(), dpi=150)
        plt.show()
        plt.close()

    def getPixelRank(self, target=None, segmap=None):
#TODO global???
        global a
        global pixels
        global location
        if not target:
            target = self.bestCandidate
        if not segmap:
            segmap = self.segmap
        a = np.where(segmap == target+1)
        pixels = []
        for k in range(len(a[0])):
            pixels.append((a[0][k], a[1][k])) #list of tuples of coords
        if not  in pixels:
            return 0
        
        def sortkey(x):
            a, b = x
            return self.swappedData[a][b]        
        pixels.sort(key = sortkey)
        location = pixels.index((int(self.event['x']), int(self.event['y'])))
        return float(location)/float(len(pixels))
    
    
    def errorProtocol(self, e, target=None):
            curFile = self.idNumString + '.' + str(self.filterNum)
            errorString = str(namecount) + '\n' + str(e) + " " + curFile + '\n'
            print(errorString)
            with open(ERRORFILE, 'a+') as errorfile:
                errorfile.write(errorString)

            chosen = target if target else self.bestCandidate
    
            if PLOT_ERR:
                self.plot(myVmin = 0, myVmax = 1000, target=chosen)
                self.plot()
            if ONLY_FLAG_ERRORS or e=='far' or e=='unlikely':
                return
            else:
                print("raising")
                raise


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


''' load event redshifts' dictionary '''
global zdict
zdict = {}
zfile = open('new_ps1z.dat', 'r')
zfile.readline() #get rid of heade

for line in zfile:
    parts = line.split()
    eventId = parts[0][3:]
    redshift = float(parts[1])
    zdict[eventId] = redshift

zfile.close()
#not in zdict, looked up from ps1 website 
#http://telescopes.rc.fas.harvard.edu/ps1/2014/alerts/
zdict['100014'] = 0.357
zdict['300220'] = 0.094
zdict['380108'] = 0.159



'''load event location dictionary'''
db = pd.read_table('alertstable_v3',sep=None,index_col = False,
               engine='python')
db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False,
                engine='python')
db = db.append(db2,ignore_index=True)



'''load prerun sdss queries'''
global fullSdssTable
fullSdssTable = Table.read('sdss_queries.dat', format='ascii')
fullSdssTable = fullSdssTable.group_by('idnum')


'''load real host data'''
global hostsData
# not sure if I should be able to json.load this, but I can't and this works
with open('hosts.dat', 'r') as hostsfile:
        hostsData = hostsfile.read()
# convert to dictionary from dictionary string
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
            if self.images[x].objCoords[self.images[x].bestCandidate].separation(event['coords']) < MINDIST:
                good_images.append(x)
            elif photozMatched:
                good_photozs.append(x)
            else:
                BAD_IMAGES.append(x)

        
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
#TODO pick in chance coincidence
            self.filter_to_use = good_photozs[0]
            goodCoords = self.images[self.filter_to_use].\
                objCoords[self.images[self.filter_to_use].bestCandidate]
            #use object at that location in all images
            for x in range(3,7):
                self.images[x].correct_bestCandidate_to(goodCoords, self.filter_to_use)
        else:
            # find minimum chance coincidence object
            first = self.images[BAD_IMAGES[0]]
            minChanceCoincidence = first.chanceCoincidence[first.bestCandidate]
            minOwner = BAD_IMAGES[0]
            for num in BAD_IMAGES:
                im = self.images[num]
                bestChance = im.chanceCoincidence[im.bestCandidate]
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
                #GET DEFAULT DATA AND RETURN
                self.images[self.filter_to_use].errorProtocol(
                        "Best.%s.NO_GOOD_CANIDIDATE" % minChanceCoincidence)
                global BAD_COUNT
                BAD_COUNT += 1
                all_sn_data = {}
                for x in range(3,7):
                    all_sn_data.update(self.images[x].getDefaultData())
                    
                self.used_default=False
                #get non filter dependent data
                all_sn_data.update(self.getSnFinalData(None))
                print('ID: %s' % self.idNum)
                return all_sn_data
            
            
            
        self.used_default = True
        all_sn_data = {}
        chosen_im = self.images[self.filter_to_use]
        # get final data from the known good filter first
        best_data = chosen_im.getFinalData()
        all_sn_data.update(best_data)
        for i in range(3,7):
            if i == self.filter_to_use: 
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
    all_redshifts = {}
    all_kronMags = {}
    all_kronRads = {}
    for snType in TYPES:
        all_redshifts[snType] = []
        all_kronMags[snType] = []
        all_kronRads[snType] = []

    all_all_data = []
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
        all_all_data.append(s.run())
        
        #TODO column order
        #print(all_all_data)
        df = pd.DataFrame(all_all_data, columns=COLUMNS)
        if WRITE_CSV:
            df.to_csv(DESTDIR + WRITE_CSV)

def main():
#TODO remove
    import time
    start = time.time()
    #figure out which files to use based on value specified at top
    if FILES == 'all' or FILES =='range' or FILES =='new random':
        #CHANGED now only looks for .3 files, assuming the rest are there
        filenames = sorted(glob.glob(SOURCEDIR + '/ps1hosts/psc*.3.fits'))
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
    print(end - start)
    print("BAD COUNT: %s" % BAD_COUNT)

if __name__ == "__main__":
     main()
     
