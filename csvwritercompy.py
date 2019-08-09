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
import json
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from matplotlib.patches import Ellipse, RegularPolygon
from astropy.io import fits
from astropy.wcs import WCS

from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table, vstack
from sdssTableGroupIndex import sdssTableGroupIndex

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
#TODO does it make sense to do limiting mag + 1? same in default final props
#TODO make note: to change a per-filter property, change in 3 places
#TODO make each property into an object?

FILTER_M0s = (None, None, None, 22.918, 22.822, 23.652, 23.540) 
# FILTER_M0s[i] is the average magnitude zero point in filter i, use by default if no stars for calibration
# Also use if star calibration yields outlier (> 1 away from average value)

#USAGE FLAGS:
WRITE_CSV = "\galaxiesdata.csv" # filename to write to or None
MAG_TEST_ALL = False
MAG_TEST_STDEV = False
PLOT_REDSHIFTS = False

PRINT_DATA = False

CHECK_DISTANCE = 5 #print all files with most likely host farther than this arcsecs
PLOT_ALL = False
PLOT_ERR =  True #plots only files that give errors or low probability
PLOT_DIR = os.getcwd() + '/plots' # where to put plot images
ONLY_FLAG_ERRORS = True # catch errors, print filename, move on
FILES = 'specified' #options are 'all', 'preset random', 'new random', 'range', 'specified', 'nonsquare

#TODO delete
SPECIFIED = []
to_check = [160103, 180313, 590123, 50296, 90034, 50601]
for f in to_check:
    SPECIFIED.extend(glob.glob((SOURCEDIR + '/ps1hosts/psc*%i*.[3-6].fits' % f)))

SPECIFIED = [SOURCEDIR + '/ps1hosts/psc480552.6.fits']
RANGE = (0,40)
m0collector = [None, None, None, [], [], [], []]
BAD_COUNT = 0
'''make header'''
COLUMNS =['ID', 'hostRa', 'hostDec', 'offby', 'hectoZ', 'redshift_dif']

perImageHeaders = ['KronRad (kpc)', 'separation (kpc)', 'area (kpc^2)', 'sep/area (kpc)',
                   'x', 'y','KronMag', 'Abs. Mag', 'Angle',
                   'Ellipticity', 'RA',  'DEC', 
                   'Discrepency (arcsecs)', 'pixelRank', 'chance coincidence',
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

class Image:    
    def __init__(self, idNumString, filterNum, event):
#TODO get all of this work out of the init? Make a run method?
        
#TODO unneccessary or redundant?
        self.idNumString = idNumString
        self.idNum = int(idNumString)
        self.filterNum = filterNum
        self.event = event
        self.eventz = float(zdict[self.idNumString])
        
    def run(self):
        filename = FILENAME_PREFIX + "%s.%s.fits" % (self.idNumString, self.filterNum)
        w = WCS(filename)
        
        
        image_file = fits.open(filename)
        
        # to get image dimensions in wcs:
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
        
        
        def recursiveExtraction(attemptedThresh):
            self.objects, self.segmap = sep.extract(self.swappedData, attemptedThresh,
                                          err=self.bkg.globalrms, #var=noise_data2 ,
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
                # remove any 'nan' kronrads and blacklist the objects
                #check if ki is nan:
                if np.isnan(self.kronrad[i]):
                    #if an object near event (likely host) fails,
                    # redo getting objects using a higher thresh
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

        ''' '''
        
        # start with default threshold, if it fails then use a higher one
        recursiveExtraction(THRESHOLD)
        
        #remove any nans, replace with cell saturation
        #TODO better solution?
        sat = image_file[0].header['HIERARCH CELL.SATURATION']
        self.swappedData[np.isnan(self.swappedData)] = sat
        
        flux = []
        for i in range(len(self.kronrad)):
            try:
                thisflux, _fluxerr, _flag = sep.sum_ellipse(self.swappedData, self.objects['x'][i],
                                        self.objects['y'][i], self.objects['a'][i],
                                        self.objects['b'][i], self.objects['theta'][i],
                                        2.5*self.kronrad[i], subpix=1)
                flux.append(thisflux)
            except:
                flux.append(0)
                self.errorProtocol("Warning: failed flux calculation on ")
        flux = np.array(flux)

        if self.idNumString not in sdssTableGroupIndex:
            sdssTable = None
        else:
            sdssTable = fullSdssTable.groups[sdssTableGroupIndex[self.idNumString]]

#TODO the last argument is 1 b/c first pixel of FITS should be 1,1 not 0,0?
#TODO icrs?
        self.ra, self.dec = w.all_pix2world(self.objects['x'], self.objects['y'], 1)
        self.objCoords = SkyCoord(self.ra*u.deg, self.dec*u.deg, frame='icrs') # coords of all objs

        # for collecting mags and fluxes to calculate the zero for this file
        colRealMags = []
        colFluxes = []
        self.photozs = [None]*len(self.objects)
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

                        # if there's a valid photz, store for future comparison
                        if sdssTable['type'][j] == 3 \
                           and sdssTable['z'][j] \
                           and not np.isnan(sdssTable['z'][j]) \
                           and sdssTable['z'][j] != -9999.0:
                               self.photozs[i] = sdssTable['z'][j]

                        # if its a star, and its not right on the event, blacklist
                        if sdssTable['type'][j] == 6 \
                            and not abs(sdssTable['ra'][j] - self.event['ra'])*u.deg < MINDIST \
                            and not abs(sdssTable['dec'][j] - self.event['dec'])*u.deg < MINDIST:
                            self.blacklist.add(i)

                        #if its a star brighter than 21, use its magnitude in zero point calculation
                        if sdssTable['type'][j] == 6 and sdssTable[FILTERS[self.filterNum]][j] < 21. \
                            and self.objects['x'][i] - 2.5 * self.kronrad[i] >= 0 \
                            and self.objects['x'][i] + 2.5 * self.kronrad[i] <= self.maxX \
                            and self.objects['y'][i] - 2.5 * self.kronrad[i] >= 0 \
                            and self.objects['y'][i] + 2.5 * self.kronrad[i] <= self.maxY:
                                colRealMags.append(sdssTable[FILTERS[self.filterNum]][j])
                                colFluxes.append(flux[i])
#TODO check that star is completely in frame before use
                        break

        colFluxes = np.array(colFluxes)

        self.exposure_time = float(image_file[0].header['EXPTIME'])
        m_0s = colRealMags + 2.5 * np.log10(colFluxes/float(self.exposure_time))
        m_0s = m_0s[~np.isnan(m_0s)] #remove nans
        if m_0s.any():
            self.m_0 = np.median(m_0s)
            global m0collector
            m0collector[self.filterNum].append(self.m_0)

        #no SDSS stars to calculate zero point
        else:
                self.m_0 = FILTER_M0s[self.filterNum]
                colFluxes = np.array([])
                colRealMags = np.array([])

        # m_0 is outlier, use default
        if abs(self.m_0 - FILTER_M0s[self.filterNum]) > 1:
                self.m_0 = FILTER_M0s[self.filterNum]


        colMyMags = -2.5 * np.log10(colFluxes/self.exposure_time) + self.m_0
        self.magnitude = -2.5 * np.log10(flux/self.exposure_time) + self.m_0
        self.magnitude[np.where(np.isnan(self.magnitude))] = LOWEST_MAG

#TODO combine?
        if MAG_TEST_ALL:
            all_myMags[self.filterNum].extend(colMyMags[:])
            all_realMags[self.filterNum].extend(colRealMags[:])
#TODO remove
            for k in range(len(colMyMags)):
                if abs(colMyMags[k] - colRealMags[k]) > 2:
                    self.errorProtocol("outlier: ")

        '''Chance Coincidence Calculation'''
        #size is half light radius
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
        # exclude any blacklisted
        for i in self.blacklist:
            self.chanceCoincidence[i] = 1
        
        self.ellipticity = 1 - (self.objects['b']/self.objects['a'])

        #remove mpc for dimensionless math and convert to kpc
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
   
    def attemptBestCandidate(self):
        self.bestCandidate = np.argmin(self.chanceCoincidence)
        # Correct bestCandidate choice for closer redshift if similar prob.
#TODO change to absolute separation
#TODO use ovals!
        photoz_matched = False
        if self.separation[self.bestCandidate] > CHECK_DISTANCE:
            if self.chanceCoincidence[self.bestCandidate] < 0.02:
                self.errorProtocol("warning, checking for matching photozs for far but likely candidate")
             
            #check for an object with photoz within 10% of event redshift
            for k in range(len(self.photozs)):
                if self.photozs[k] and abs(self.photozs[k] - self.eventz)/self.eventz < 0.1:
                    if photoz_matched:
                        self.errorProtocol("multiple matching photozs")
                        if self.chanceCoincidence[k] < self.chanceCoincidence[self.bestCandidate]:
                               self.bestCandidate = k
                        photoz_matched = True

        #check for objects cut off
        if self.objects['x'][self.bestCandidate] + 2.5* self.kronradKpc[self.bestCandidate] > self.maxX or \
            self.objects['x'][self.bestCandidate] - 2.5* self.kronradKpc[self.bestCandidate] < 0 or \
            self.objects['y'][self.bestCandidate] + 2.5* self.kronradKpc[self.bestCandidate] > self.maxY or \
            self.objects['y'][self.bestCandidate] - 2.5* self.kronradKpc[self.bestCandidate] < 0:
                self.errorProtocol("cut off object")                        
                        
#TODO green circle all in plot
        return photoz_matched
    
    def correct_bestCandidate_to(self, goodcoords, used_filter):
        if used_filter.filterNum == self.filterNum:
            return
        self.bestCandidate = None
        for k in range(len(self.objCoords)):
            if self.objCoords[k].separation(goodcoords) < MINDIST:
                self.bestCandidate = k
#TODO return true or false, check
    
    def modifyFinalDataFrom(self, data, newFluxParams, used_bestCandidate, used_segmap):
        # no objects in correct location detected in bad image, 
        #use data from good image but update magnitude
        newFlux, _fluxerr, _flag = sep.sum_ellipse(*newFluxParams, subpix=1)
        newMagnitude = -2.5 * np.log10(self.exposure_time) + self.m_0
#TODO just save the old pixel rank
        newPixelRank = self.getPixelRank(target=used_bestCandidate, segmap=used_segmap)
        newAbsMag = newMagnitude - 5*np.log10(self.dL) - 10
        
        f = self.filternum 
        old_filternum = data.keys()[0][-1]
        oldMagnitude = data['KronMag_%s' % old_filternum]
        new_data = {}
        for (k, v) in data.iteritems():
            new_key = k[:-1] + str(self.filterNum)
            new_data[new_key] = v
        new_data['KronMag_%s' % f] = newMagnitude
        new_data['Abs. Mag_%s' % f] = newAbsMag
        new_data['pixelRank_%s' % f] = newPixelRank
        new_data['KronRad (kpc)_%s' % f] = newMagnitude
        oldChanceCoincidence = data['chanceCoincidence_%s' % old_filternum]
#TODO check length, make sure this works
        #adjust chanceCoincidence for change in magnitude
        new_data['chanceCoincidence_%s' % f] = 1 - (1 - oldChanceCoincidence)**(10**(0.33*(newMagnitude - oldMagnitude)))
#TODO make sure above math is correct
        new_data['host_found_%s' % f] = 0
    
    def getFinalData(self):    
        '''Get "real" host location and redshifts'''
        finalDict = {}
        f = self.filternum               
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
    return finalProperties

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
            #TODO note that redshift_offset of -1 indicates no photoz, positive is difference
            defaultFinalProperties['host_found_%s' % f] = 0
    return defaultFinalProperties
   
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
        if not (int(self.event['x']), int(self.event['y'])) in pixels:
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
        
    def getSnFinalData(self):
        finalDict = {}
        hostRa = hostsData[self.idNumString]['host_ra']
        hostDec = hostsData[self.idNumString]['host_dec']
        try:
            hostCoords = SkyCoord(hostRa, hostDec, unit=(u.hourangle, u.deg))
            #convert to decimal deg:
            finalDict['hostRa'] = hostCoords.ra.deg
            finalDict['hostDec'] = hostCoords.dec.deg
            finalDict['offby'] = hostCoords.separation(self.objCoords[self.bestCandidate]).arcsec\
                    if hostRa else None
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
       finalDict['redshift_dif'] = -1
        
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
            goodCoords = self.images[self.filter_to_use].objCoords[self.images[self.filter_to_use].bestCandidate]
#TODO do I use object at event coords or good coords?
            for x in good_photozs.extend(BAD_IMAGES):
                self.images[x].correct_bestCandidate_to(goodCoords, self.filter_to_use)
        elif good_photozs:
#TODO pick in chance coincidence
            self.filter_to_use = good_photozs[0]
            goodCoords = self.images[self.filter_to_use].objCoords[self.images[self.filter_to_use].bestCandidate]
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
            self.filter_to_use=num
            if minChanceCoincidence <= 0.02:
                self.images[self.filter_to_use].errorProtocol("USING BEST CHANCE: %s" % minChanceCoincidence)
                goodCoords = self.images[self.filter_to_use].objCoords[self.images[self.filter_to_use].bestCandidate]
                for i in range(3,7):
                    self.images[i].correct_bestCandidate_to(goodCoords, self.filter_to_use)
            else:
                self.images[self.filter_to_use].errorProtocol("Best.%s.NO_GOOD_CANIDIDATE" % minChanceCoincidence)
                global BAD_COUNT
                BAD_COUNT += 1
                for x in range(3,7):
                    row.update(self.images[x].getDefaultData())
                    return

        all_sn_data = {}
        # get final data from the known good filter first
        best_data = self.images[self.filter_to_use].getFinalData()
        all_sn_data.update(best_data)
        for i in range(3,7):
            if i == self.filter_to_use: #alread done above
                pass
            elif self.images[i].bestCandidate == None:
#TODO check correctnesss
                usedbestCandidate = self.images[self.filter_to_use].bestCandidate
                newFluxParams =  (self.images[i].swappedData,
                   self.images[self.filter_to_use].objects['x'][usedbestCandidate],
                   self.images[self.filter_to_use].objects['y'][usedbestCandidate],
                   self.images[self.filter_to_use].objects['a'][usedbestCandidate],
                   self.images[self.filter_to_use].objects['b'][usedbestCandidate],
                   self.images[self.filter_to_use].objects['theta'][usedbestCandidate], 
                   2.5*self.images[self.filter_to_use].kronradKpc[usedbestCandidate])
                all_sn_data.update(self.images[i].modifyFinalDataFrom(bestData, newFluxParams))
            else:
                all_sn_data.update(self.images[i].getFinalData())
#TODO make flattening better                
        #return (cache[3].extend(cache[4].extend(cache[5].extend(cache[6]))))
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

    if WRITE_CSV:
        #create destination directory if it does not exist
        if not os.path.isdir(DESTDIR):
            os.mkdir(DESTDIR)
        destfile = open(DESTDIR + WRITE_CSV, "w+")

    if PRINT_DATA:
        print(HEADER)
    all_all_data = []
    for filename in filenames: 
        # to extract transient and filter number from filename of the form
        # */psc190369.6.fits
        dotSplit = filename.split('.')
        # transient identifier is everything after the last 'c' but before
        # the second to last dot
        idNumString = dotSplit[-3].split('c')[-1] #used for getting hectospec data
        # filter number is between second to last dot and last dot
        s = Supernova(idNumString)
        all_all_data.append(s.run())
        #TODO column order
        df = pd.DataFrame(all_all_data, columns=COLUMNS)
        if WRITE_CSV:
            df.to_csv(WRITE_CSV)
 
    if WRITE_CSV:
        destfile.close()
    if PRINT_DATA:
        print(all_sn_data)

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
     
