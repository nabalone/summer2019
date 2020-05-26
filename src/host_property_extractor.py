# -*- coding: utf-8 -*-

#TODO* check that all constants are at top, not hardcoded in

'''
fits files must be named pscxxxxxx.f.fits
and FILTERS must be such that filter number f corresponds with filter filters[f]
alertstable_vs and alertstablevs.lasthalf must be in directory for sn locations
sdss_queries.dat

'''

import numpy as np
import pandas as pd
import os
import glob
import random
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
import json

#attempt to supress fits warnings about date changed, etc.
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

PROJ_HOME = os.environ['DATA_SRCDIR']

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--mask', action='store_true', 
                    help='only generate mask and date files for cnn')
parser.add_argument('--use_prev', action='store_true', 
                    help='This is not the first time script has been run; \
                        use mag zero points and property averages from last run')
args, _remaining = parser.parse_known_args() #args = parser.parse_args()
if args.mask:
    print("mask making only")
if args.use_prev:
    print("Using m_0s and averages from a prev. run")

SOURCEDIR = PROJ_HOME + '/src/all_fits' #fits images location
DICTDIR = PROJ_HOME + '/src' #all other data files location
OUTPUT_DIR = PROJ_HOME + '/src/outputs'
DESTDIR = os.getcwd()

ERRORFILE = OUTPUT_DIR + '/errorfile.txt' #really just a log file of the run
WRITE_CSV = OUTPUT_DIR + "/galaxiesdata.csv" # filename to write to or None
PLOT_DIR = OUTPUT_DIR + '/plots' # where to put plot images
PLOT_DIR2 = OUTPUT_DIR + '/zsorted_plots' #redshift-sorted plots for paper
# clear leftover plots from previous runs
old_plots = glob.glob(PLOT_DIR + '/*') + glob.glob(PLOT_DIR2 + '/*')
for old_plot in old_plots:
    os.remove(old_plot)
if not os.path.isdir(PLOT_DIR): 
    os.mkdir(PLOT_DIR)
if not os.path.isdir(PLOT_DIR2): 
    os.mkdir(PLOT_DIR2)
FILENAME_PREFIX = SOURCEDIR + "/psc" #everything before the sn number

FILLER_VAL = None
THRESHOLD = 3 #sigma of detection
MAXTHRESH = 30 # Do not raise threshhold beyond this point. somewhat arbitrary.
PSF = 4 #the FWHM
MINAREA = 3.14 * (PSF/2)**2 #area of a circle with psf diameter
DEBLEND_CONT = 0.01 # for sep.extract. 1.0 turns off deblending, 0.005 is default
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as
FILTERS = [None, None, None, 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z']
#in fits file naming psc12345.x.fits, filter numbers x of 3,4,5,6 correspond to
#filters g,r,i,z respectively
LIKELIHOOD_THRESH = 0.2 # only choose hosts with chance coincidence below this 
TYPES = {'SNIIn', 'SNIa', 'SNII', 'SNIbc', 'SLSNe'}

#has been replaced by two_sigma_mag
#Magnitude used for all effectively undetectable objects
#LOWEST_MAG = 26 #the limiting mag is 25

#MAG_TEST_ALL = False

#USAGE FLAGS that I never put into the argparser:
PLOT_ALL = True
PLOT_ERR =  True #plots files when something eventful happens e.g. errors, low probabilities
ONLY_FLAG_ERRORS = True # catch errors, log and move on

#TODO restore
FILES = 'range' #options are 'all', 'new random', 'range', 
#'specified', 'nonsquare'. Use 'all' for normal usage, the rest are for debugging.
#if 'specified', will run on the subset of files specified in SPECIFIED. 
#if 'range', will run on subset from RANGE[0] file to the RANGE[1] file
#e.g.:
SPECIFIED = [SOURCEDIR + '/psc020121.3.fits', 
             SOURCEDIR + '/psc020121.4.fits',
             SOURCEDIR + '/psc020121.5.fits',
             SOURCEDIR + '/psc020121.6.fits']
RANGE = (3,8)


m0collector = [None, None, None, [], [], [], []] #saved for use in future runs
BAD_COUNT = 0 #number of sn with no good host detected


# IMPORTANT NOTE: to change what properties are in the excel file, 
# make the change in the header, in the dictionaries 
# finalDict in Image.getFinalData and new_data in Image.modifyFinalDataFrom
# or Supernova.getSnFinalData
# and change defaults in Image.getDefaultData and the collection of previous
# average data i.e. setting AVRG_* (e.g. AVRG_RADII) if appropriate. 
'''make header'''
COLUMNS =['ID', 'Redshift', 'comparison_hectoZ', 
          'Event RA', 'Event DEC', 'Host RA', 'Host DEC', 
          'comparison_hostRa', 'comparison_hostDec', 
          'comparison_distance_dif (arcsec)']

perImageHeaders = ['KronRad (kpc)', 'separation (kpc)', 'area (kpc^2)', 
                   'sep/sqrt(area) (kpc)',
                   'x', 'y','KronMag', 'Abs. Mag', 'Angle',
                   'Ellipticity', #'RA',  'DEC', 
                   #'Discrepency (arcsecs)', 
                   'pixelRank', #'redshift', 
                   'chanceCoincidence',
                   'host_found']
for i in range(3,7):
    for val in perImageHeaders:
        COLUMNS.append(val + '_' + str(i))

# clear any previous logs from errorfile
if os.path.exists(ERRORFILE):
  os.remove(ERRORFILE)
  
  
#all_myMags = []
#all_realMags = []

#for checking how many sn are actually inside their host
# just out of curiosity
#INSIDE = {}
#OUTSIDE = {}
#INSIDE['SNIa'] = set()
#INSIDE['SNIbc'] = set()
#INSIDE['SNII'] = set()
#INSIDE['SNIIn'] = set()
#INSIDE['SLSNe'] = set()
#
#OUTSIDE['SNIa'] = set()
#OUTSIDE['SNIbc'] = set()
#OUTSIDE['SNII'] = set()
#OUTSIDE['SNIIn'] = set()
#OUTSIDE['SLSNe'] = set()


# generates sequential numbers for uniquely naming plots files
namecount = 0
def namecountgen():
    global namecount
    namecount += 1
    return namecount


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
    
        
    '''load "real" host galaxy data'''
    # NOT calculated by this script. For testing, comparison, to check our work
    # not used for actual research/results
    global comparisonHostsData
    with open(DICTDIR + '/hosts.dat', 'r') as hostsfile:
             comparisonHostsData = hostsfile.read()
# not sure if I should be able to json.load this, but it fails and this works
    comparisonHostsData = ast.literal_eval(comparisonHostsData)



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
        image_file = fits.open(filename)
        
        # to get image corners
        self.maxX = image_file[0].header['NAXIS1']
        self.maxY = image_file[0].header['NAXIS2']
        
        ''' extract objects '''
        image_data = image_file[0].data

        #fix byte order
        self.swappedData = image_data.byteswap(True).newbyteorder()
        # subtracting out background
        self.bkg = sep.Background(self.swappedData)
        self.swappedData = self.swappedData - self.bkg      
        
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
            # Can also be used for detecting if event location is in an object
            # if the code involving the variables INSIDE and OUTSIDE is uncommented
            self.objects, self.segmap = sep.extract(self.swappedData, attemptedThresh,
                                          err=self.bkg.rms(), 
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
                # An objcect may have a 'nan' kronrad for a few reasons. 
                # Usually because [I forgot reason, check and justify ignoring]
#TODO
                # Or could be an oversaturated object 
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
                    # which will not be used anyway
                    self.blacklist.add(i)
                    self.kronrad[i] = 0.1

        # start with default threshold
        recursiveExtraction(THRESHOLD)
        
        # remove any nan pixels in image, which are due to saturation
        # replace with cell saturation value for closest magnitude accuracy
        # this must be done after recursiveExtraction to make sure we are 
        # separating host from oversaturated objects
        #TODO is there a better solution than using cell saturation level?
        sat = image_file[0].header['HIERARCH CELL.SATURATION']
        self.swappedData[np.isnan(self.swappedData)] = sat
        
        flux = []
        for i in range(len(self.kronrad)):
            try:
                thisflux, _fluxerr, _flag = self.sum_ellipse_extrapolated(i)
                flux.append(thisflux)
            except:
                # probably majority of object was outside frame or object was
                # too narrow; unlikely to be host, ignore object
                self.blacklist.add(i)
                flux.append(0)
                #deprecated logging:
                #self.bestCandidate = i
                #self.errorProtocol("Warning: failed flux calculation on ")
        flux = np.array(flux)

        # fullSdssTable should contain pre-saved querries of event locations
        sdssTable = fullSdssTable[fullSdssTable['idnum']==self.idNum]
            
        # get WCS coords of all detected objects
        w = WCS(filename)
        self.ra, self.dec = w.all_pix2world(self.objects['x'], self.objects['y'], 0)
        self.objCoords = SkyCoord(self.ra*u.deg, self.dec*u.deg, frame='icrs') # coords of all objs

        # for collecting mags and fluxes to calculate the zero point for this file
        colRealMags = []
        colFluxes = []
        
#TODO* check that magnitudes are correct
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

#TODO check if two_sigma value is reasonable
        # the flux of a 2 sigma detection in the least-noisy part        
        # the lowest possible two sigma value in the image times the PSF area
        two_sigma_flux = 2 * np.min(self.bkg.rms()) * (np.pi * PSF/2)
        self.two_sigma_mag = -2.5 * \
                np.log10(two_sigma_flux/self.exposure_time) + self.m_0
        
        # calculate all magnitudes using the new zero point
        self.magnitude = -2.5 * np.log10(flux/self.exposure_time) + self.m_0
        
        # replace all nan magnitudes with two_sigma detection magnitude
        self.magnitude[np.where(np.isnan(self.magnitude))] = self.two_sigma_mag
            

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
        self.absLimMag = self.two_sigma_mag - 5*np.log10(self.dL) - 10 + 2.5 * np.log10(1.+self.eventz)
        image_file.close()
        
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

        if not self.isEventIn(self.bestCandidate):
            if self.chanceCoincidence[self.bestCandidate] < 0.02:
                self.errorProtocol("warning, checking for matching photozs for far but likely candidate")
             
            # check for objects with photoz within 10% of event redshift
            # choose matching redshift object with lowest chance coincidence
            for k in range(len(self.photozs)):
                if self.photozs[k] and abs(self.photozs[k] - self.eventz)/self.eventz < 0.1:
                    if not photoz_matched: #first time we find matching photoz
                        #use instead
                        self.bestCandidate = k
                        photoz_matched = True
                    else: #this is 2nd+ matching photoz galaxy we have found
                        self.errorProtocol("multiple matching photozs")
                        #use instead only if lower chanceCoincidence 
                        if self.chanceCoincidence[k] < self.chanceCoincidence[self.bestCandidate]:
                               self.bestCandidate = k

                  
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
#TODO* make sure this is correct: old pixel rank used
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
#TODO* check length, make sure this works
        #adjust chanceCoincidence for the new magnitude
        new_data['chanceCoincidence_%s' % f] = 1 - (1 - oldChanceCoincidence)**(10**(0.33*(newMagnitude - oldMagnitude)))
#TODO* make sure above math is correct
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
        #finalDict['RA_%s' % f] = self.ra[bestCandidate]
        #finalDict['DEC_%s' % f] = self.dec[bestCandidate]
#TODO was Discrepency supposed to be something?      
        #finalDict['Discrepency (arcsecs)_%s' % f] = self.dec[bestCandidate]
        
        
        finalDict['pixelRank_%s' % f] = self.getPixelRank()
        finalDict['chanceCoincidence_%s' % f] = self.chanceCoincidence[bestCandidate]
        finalDict['host_found_%s' % f] = 1
        #finalDict['redshift_%s' % f] = self.eventz
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
            defaultFinalProperties['KronMag_%s' % f] = self.two_sigma_mag
            
#TODO**            FIXXX to 2sigma
            
            defaultFinalProperties['Abs. Mag_%s' % f] = self.absLimMag
            if args.use_prev:
                defaultFinalProperties['Angle_%s' % f] = AVRG_ANGLES[self.filterNum]
                defaultFinalProperties['Ellipticity_%s' % f] =  AVRG_ELLIPTICITIES[self.filterNum]
            else:
                defaultFinalProperties['Angle_%s' % f] = 0
                defaultFinalProperties['Ellipticity_%s' % f] =  0.3
            #defaultFinalProperties['RA_%s' % f] = self.event['ra']
            #defaultFinalProperties['DEC_%s' % f] = self.event['dec']
            #defaultFinalProperties['Discrepency (arcsecs)_%s' % f] = None
            if args.use_prev:                
                defaultFinalProperties['pixelRank_%s' % f] = AVRG_PIXELRANKS[self.filterNum]
            else:
                defaultFinalProperties['pixelRank_%s' % f] = 0.5
            defaultFinalProperties['chanceCoincidence_%s' % f] = 1
            defaultFinalProperties['host_found_%s' % f] = 0
#            if args.use_prev:
#                defaultFinalProperties['redshift_%s' % f] = AVRG_REDSHIFTS[self.filterNum]
#            else:
#                defaultFinalProperties['redshift_%s' % f] = 0.2
        return defaultFinalProperties


#TODO remove, unused
    # returns |Redshift_SN - Photoz_galaxy|
    # or -1 if no SDSS photoz best candidate galaxy
    # we don't have enough sig.figs to compare with photozerr
#    def getPhotozDif(self):
#        photoz = self.photozs[self.bestCandidate]
#        eventz = self.eventz
#        
#        #PHOTOZ failed or no match
#        if photoz == None: 
#            return -1
#        else:
#            return abs(eventz - photoz)
        

    # saves and displays image with SDSS identified stars and 
    # blacklisted (disqualified) objects circled in blue,
    # best host candidate circled in green, all other detected objects circled
    # in red, event location marked with purple triangle. Saves to PLOT_DIR
    def plot(self, myVmin=None, myVmax=None, target=None, my_title=''):
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
            e.set_facecolor('none')
            if i in green:
                e.set_edgecolor('green')
            elif i in self.blacklist:
                e.set_edgecolor('blue')
            else:
                e.set_edgecolor('red')
            ax.add_artist(e)

        if my_title=='final':
        # the plots going into PLOT_DIR are for debugging; only want 1 filter
        # plots going into PLOT_DIR2 are for the paper, so we want all filters
            if self.filterNum==5:
                    plt.title(my_title)
                    plt.savefig(PLOT_DIR + "/galaxyimage" + self.idNumString + '_' \
                    + str(namecountgen()) + ".png", dpi=150)
            plt.title(' ')
            filterletter={3:'g', 4:'r', 5:'i',6:'z'}
            plt.savefig(PLOT_DIR2 + "/z%s_sn%s_type%s_%s.png" % (self.eventz, \
                        self.idNumString, typeDict[self.idNumString], \
                        filterletter[self.filterNum]))
        else:
            if self.filterNum==5:
                plt.title(my_title)
                plt.savefig(PLOT_DIR + "/galaxyimage" + self.idNumString + '_' \
                        + str(namecountgen()) + ".png", dpi=150)
        plt.show()
        plt.close()

    def getPixelRank(self, target=None, segmap=None):
        # NOTE SEP'S SEGMAP INDEXES Y-FIRST; 
        # OBJECT N IS CENTERED AROUND segmap[objects['y'][N]][objects['x'][N]]
        if not target:
            target = self.bestCandidate
        if not segmap:
            segmap = self.segmap
        #"All pixels belonging to the i-th object (e.g., objects[i]) have value i+1"
        #https://sep.readthedocs.io/en/v1.0.x/api/sep.extract.html
        a = np.where(segmap == target+1) 
        pixels = []
        for k in range(len(a[0])):
            pixels.append((a[0][k], a[1][k])) #list of tuples of coords
#            global INSIDE
#            INSIDE[typeDict[self.idNumString]].add(self.idNumString)
            
        #If sn occurs outside host, pixel rank is automatically 0
        if not (int(self.event['y']), int(self.event['x'])) in pixels:
#            global OUTSIDE
#            OUTSIDE[typeDict[self.idNumString]].add(self.idNumString)
            return 0
        
        #sort pixels and see what fraction of them are dimmer than the sn location pixel
        def sortkey(k):
            a, b = k
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
                self.plot(myVmin = 0, myVmax = 3000, target=chosen, my_title=my_title)
                #self.plot()
            if ONLY_FLAG_ERRORS or e=='far' or e=='unlikely':
                return
            else:
                print("raising")
                raise
                
                

class Supernova:
    
    def __init__(self, idNumString):
        self.idNumString = idNumString
        self.idNum = int(idNumString)
        
        
    # comparisonHostsData and anything starting with 'comparison_' is for comparison with 
    # data from the hosts.dat file, NOT calculated by this script
    # and will be used for testing, comparison, to check our work
    # not used for actual research/results
    def getSnFinalData(self, chosen_loc): #chosen_loc is location of our chosen host
        best_image = self.images[self.filter_to_use]
        finalDict = {'ID':self.idNum}
        finalDict['Redshift'] = best_image.eventz
        finalDict['Event RA'] = best_image.event['ra']
        finalDict['Event DEC'] = best_image.event['dec']
        finalDict['Host RA'] = best_image.ra[best_image.bestCandidate]
        finalDict['Host DEC'] = best_image.dec[best_image.bestCandidate]
        
        comparison_hostRa =  comparisonHostsData[self.idNumString]['host_ra'] #'real'
        comparison_hostDec = comparisonHostsData[self.idNumString]['host_dec'] #'real'
        try:
            #convert to decimal deg:
            hostCoords = SkyCoord(comparison_hostRa, comparison_hostDec, unit=(u.hourangle, u.deg))
            finalDict['comparison_hostRa'] = hostCoords.ra.deg #'real'
            finalDict['comparison_hostDec'] = hostCoords.dec.deg #'real'
            if comparison_hostRa and chosen_loc:
                # distance between our chosen host and theirs
                #TODO change comparison_distance_dif to kpc
                finalDict['comparison_distance_dif (arcsec)'] = hostCoords.separation(chosen_loc).arcsec
            else:
                # no 'real' host data for this sn
                finalDict['comparison_distance_dif (arcsec)'] =None
            finalDict['comparison_hectoZ'] = comparisonHostsData[self.idNumString]['redshift'] #'real'
        except ValueError:
            if len(comparison_hostRa) > 12: # hecto gave multiple host galaxies
                finalDict['comparison_hostRa'] = "multiple"
                finalDict['comparison_hostDec'] = "multiple"
                finalDict['comparison_distance_dif (arcsec)'] = None
                finalDict['comparison_hectoZ'] = "multiple"
            else:
                raise

        return finalDict
        
    def run(self):
        # find pixel coordinates of event; not all images are square
        e = db.where(db['eventID'] == self.idNum).dropna()
        eRa = e['ra'].values[0] #values gives np arrays
        eDec = e['dec'].values[0] #'hh:mm:ss.sss'
        #TODO: check coordinate frame consistency, ircs or fk5?
        
        #save a dict of data of event (sn) location
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
        good_images = [] #images in which a host surrounding sn location was found
        BAD_IMAGES = []  #images in which no good host was found
        good_photozs = [] #images in which a host was found farther but with matching redshift
        for x in range(3,7):
            self.images[x] = Image(self.idNumString, x, event)
            self.images[x].run()
            
            #computes and stores the best host candidate, then returns whether
            #it had a photoz matching the sn redshift
            photozMatched = self.images[x].attemptBestCandidate()
            
            #check if the host chosen for this filter has the sn within it
            #if so put this filter number good_images
            chosen_host = self.images[x].bestCandidate
            if chosen_host != None: 
                if self.images[x].isEventIn(chosen_host):
                    good_images.append(x)
#kronrad is highly variable, don't use here
#                separation = self.images[x].objCoords[chosen_host].separation(event['coords'])
#                kronrad = self.images[x].kronrad[chosen_host]
#                if separation < kronrad:
#                    good_images.append(x)
                    
                elif photozMatched:
                    #chosen host was found farther but had matching redshift
                    good_photozs.append(x)
                    
                else: #no good candidate 
                    BAD_IMAGES.append(x)
                
            else: #no good candidate 
                BAD_IMAGES.append(x)

        
        if args.mask:
            global masks
        
        '''choose candidates'''  
        #TODO* rewrite for clarity
        '''compare between filters, looking for the ultimate host of the 
        filters which found a host surrounding the sn location, pick the 
        ultimate host from the filter with the lowest chance coincidence. 
        (Even though they are probably all the same). For all the filters 
        who did NOT find and choose a host surrounding the SN, force them to 
        choose a host matching the ultimate host.
        
        If however no filter chose a host which was surrounding the sn, 
        look at the filters which chose a host with matching photoz to the 
        sn redshift. Of those, take the host with the lowest chance coincidence
        to be the ultimate host. Correct all other filters to match the 
        ultimate host.

        If no hosts with matching photoz were found either, take the host from 
        all the filters which had the lowest chance coincidence, provided that 
        chance coincidence was less than 0.2, and force all the other filters 
        to match it.

        However if no host with chance coincidence under 0.2 was found in any 
        of the filters, use default data
        '''
        
        def lowest_cc_of(IMAGE_LIST):
            # find minimum chance coincidence object
            if len(IMAGE_LIST) == 0:
                raise Exception("lowest_cc_of function should not be called an empty list")
            first = self.images[IMAGE_LIST[0]]
            if len(first.objects) > 0:
                minChanceCoincidence = first.chanceCoincidence[first.bestCandidate]
            else:
                # no objects were detected in image
                minChanceCoincidence = 1
            minOwner = IMAGE_LIST[0]
            for i in IMAGE_LIST:
                im = self.images[i]
                if len(im.objects) > 0:
                    i_bestChance = im.chanceCoincidence[im.bestCandidate]
                else: # no objects were detected in image
                    i_bestChance = 1
                if i_bestChance < minChanceCoincidence:
                   minChanceCoincidence = i_bestChance
                   minOwner = i
            return (minOwner, minChanceCoincidence)
        
        #filter_to_use is the filter from which host is officially chosen
        self.used_default=False #will be set to true if no decent images
        if good_images:
            self.filter_to_use, _minChanceCoincidence = lowest_cc_of(good_images)
            goodCoords = self.images[self.filter_to_use].\
                objCoords[self.images[self.filter_to_use].bestCandidate]
            #use object at that location in all images
            for x in range(3,7):
                self.images[x].correct_bestCandidate_to(goodCoords, self.filter_to_use)
            #for x in good_photozs + BAD_IMAGES:
            #    self.images[x].correct_bestCandidate_to(goodCoords, self.filter_to_use)
        elif good_photozs:
            self.filter_to_use, minChanceCoincidence = lowest_cc_of(good_photozs)
            goodCoords = self.images[self.filter_to_use].\
                objCoords[self.images[self.filter_to_use].bestCandidate]
            #use object at that location in all images
            for x in range(3,7):
                self.images[x].correct_bestCandidate_to(goodCoords, self.filter_to_use)
            self.images[self.filter_to_use].errorProtocol(
                "Using photoz matched with cc: %s" % minChanceCoincidence)
#TODO* test all combinations of good images, bad_images, no objects images            
        else:
            self.filter_to_use, minChanceCoincidence = lowest_cc_of(BAD_IMAGES)
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
                self.used_default=True
                #get non-filter-dependent data
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
            
        all_sn_data = {}
        chosen_im = self.images[self.filter_to_use]
        # get final data from the known good filter first
        best_data = chosen_im.getFinalData()
        all_sn_data.update(best_data)
        for i in range(3,7):
            if i == self.filter_to_use: 
                if args.mask: # only collecting mask
                    mask = self.images[i].segmap
                    mask = np.where(mask == self.images[i].bestCandidate + 1, 
                                    1 - self.images[i].chanceCoincidence[self.images[i].bestCandidate], 0)
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
#TODO* check correctnesss
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
        
    
        if PLOT_ALL and not self.used_default:
            for image in self.images[3:7]:
            #image = self.images[5]
                image.plot(myVmin = 0, myVmax = 3000, my_title='final') 
                
        return all_sn_data


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
        global PLOT_ALL
        global PLOT_ERR
        PLOT_ALL = False
        PLOT_ERR = False
        global ERRORFILE
        ERRORFILE = OUTPUT_DIR + '/errorfile_maskmaking.txt'
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
#TODO idNumString parsing could be done more robustly?
        
        if args.mask:
            print(filename)
            # the run notices args.mask and will save masks to the global masks variable
            s.run()
            
        else:            
            all_all_data.append(s.run())
            
            df = pd.DataFrame(all_all_data, columns=COLUMNS)
            if WRITE_CSV:
                df.to_csv(WRITE_CSV)
                
    if args.mask:
        np.savez(OUTPUT_DIR + '/all_masks', **masks)
        print('saved all masks')

            

def main():
    import time
    start = time.time()
    
    '''collect defaults from previous run if applicable'''
    #This script (host_property_extractor.py) should be run twice so that
    #the average magnitude zero points and property values from the first run can 
    #be used in the second. use without --use_prev flag the first time, and with
    # --use_prev flag the second. 
    global FILTER_M0s
    global AVRG_OFFSETS
    global AVRG_RADII 
    global AVRG_ANGLES 
    global AVRG_PIXELRANKS 
    global AVRG_ELLIPTICITIES 
    global AVRG_REDSHIFTS 
    if args.use_prev: 
        #collect average m0 from previous run for each filter, to be used as defaults
        if os.path.isfile(OUTPUT_DIR + '/m0collector.npy'):
            prev_m0s = np.load(OUTPUT_DIR + '/m0collector.npy', allow_pickle=True)
            FILTER_M0s = [0]*3
            for i in prev_m0s[3:]:
                    FILTER_M0s.append(np.median(i))
        else:
            raise Exception("Could not find m0collector.npy from previous run. \
                                Run again without --use_prev flag?")
        #collect median properties from previous run, to be used as defaults for images where 
        #no host can be detected
        #TODO I'm not sure if using means would be better than medians. 
        #There may be very bad outliers...
        if os.path.isfile(WRITE_CSV):
            prev_csv = pd.read_csv(WRITE_CSV)
            AVRG_OFFSETS = [None]*7
            AVRG_RADII = [None]*7
            AVRG_ANGLES = [None]*7
            AVRG_PIXELRANKS = [None]*7
            AVRG_ELLIPTICITIES = [None]*7
            AVRG_REDSHIFTS = [None]*7
            
            good_rows = np.array(np.where(prev_csv['host_found_3']==1))
            def get_good_med(header):
                arr = np.array(prev_csv[header])
                return np.median(arr[good_rows])
            for i in range(3,7):
                AVRG_OFFSETS[i] = get_good_med('separation (kpc)_%s'%i)
                AVRG_RADII[i] =  get_good_med('KronRad (kpc)_%s'%i)
                AVRG_ANGLES[i] =  get_good_med('Angle_%s'%i)
                AVRG_PIXELRANKS[i] =  get_good_med('pixelRank_%s'%i)
                AVRG_ELLIPTICITIES[i] =  get_good_med('Ellipticity_%s'%i)
                AVRG_REDSHIFTS[i] =  get_good_med('Redshift')
    
        else:
            raise Exception("Count not find %s from previous run. \
                                Run again without --use_prev flag?")
    else:
        print("using hardcoded default m0s")
        FILTER_M0s = (None, None, None, 22.918, 22.822, 23.652, 23.540) 
    # FILTER_M0s[i] is the average magnitude zero point in filter i, 
    #use by default if no stars for calibration
    # Also use if star calibration yields outlier (> 1 away from average value)
    
    
    
    #figure out which files to use based on value specified at top
    #'all' is for science use, the other options are just for debugging
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
    np.save(OUTPUT_DIR + "/m0collector", m0collector)

if __name__ == "__main__":
     main()
     
