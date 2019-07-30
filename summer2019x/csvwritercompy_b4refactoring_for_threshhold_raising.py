# -*- coding: utf-8 -*-

#TODO: all constants at top
#TODO: vertical line
'''
fits files must be named pscxxxxxx.f.fits
and FILTERS must be such that filter number f corresponds with filter filters[f]
alertstable_vs and alertstablevs.lasthalf must be in directory for sn locations

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
from astroquery.sdss import SDSS

SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd()
FILLER_VAL = None
THRESHOLD = 1.5
MINAREA = 5
DEBLEND_CONT = 0.01 # for sep.extract. 1.0 to turn off deblending, 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as
FILTERS = [None, None, None, 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z']
LIKELIHOOD_THRESH = 0.2
TYPES = ['SNIIn', 'SNIa', 'SNII', 'SNIbc', 'SLSNe']
TYPE_COLORS = {'SNIIn':'co', 'SNIa':'ro', 'SNII': 'bo', 'SNIbc':'go', 'SLSNe': 'mo'}
# meds used for m_0 if no sdss stars in image, or if m_0 is an outlier 
# meaning outside of upper/lower bounds
#FILTER_M0s = [None, None, None, 
#              {'upper': 23.625, 'default': 23.125, 'lower': 22.625}, # filter 3
#              {'upper': 23.4, 'default': 22.9, 'lower': 22.4}, # filter 4
#              {'upper': 23.4, 'default': 22.9, 'lower': 22.4}, # filter 5
#              {'upper': 24.125, 'default': 23.625, 'lower': 23.125}] # filter 6 

#update from first 120 files:
FILTER_M0s = [None, None, None, 
              {'upper': 23.54, 'default': 23.04, 'lower': 22.54}, # filter 3
              {'upper': 23.49, 'default': 22.99, 'lower': 22.49}, # filter 4
              {'upper': 24.29, 'default': 23.79, 'lower': 23.29}, # filter 5
              {'upper': 24.23, 'default': 23.73, 'lower': 23.23}] # filter 6 


#USAGE FLAGS:
WRITE_CSV = "\galaxiesdata.csv" # filename to write to or None
MAG_TEST_ALL = False
MAG_TEST_STDEV = False
PLOT_REDSHIFTS = False

MAGFILE = 'newmagtest.csv'
#Do not write csv an mag test simultaneously because there is only 1 writer
#TODO separate the csvwriters
#ACTUALLY< REMOVING MAG WRITER
#if MAG_TEST_ALL or MAG_TEST_STDEV:
#    WRITE_CSV = False
PRINT_DATA = False

CHECK_DISTANCE = 5 #print all files with most likely host farther than this arcsecs
PLOT_ALL = False
PLOT_ERR =  True #plots only files that give errors or low probability
PLOT_DIR = os.getcwd() + '/ploterr120' # where to put plot images
ONLY_FLAG_ERRORS = True # catch errors, print filename, move on
FILES = 'specified' #options are 'all', 'preset random', 'new random', 'range', 'specified', 'nonsquare
SPECIFIED = ['C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc000190.4.fits']
            #'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc050081.4.fits'] 
             #'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc091040.3.fits',
             #'C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc130945.5.fits']
#[SOURCEDIR + '\psc000010.3.fits', SOURCEDIR + '\psc000010.4.fits', SOURCEDIR + '\psc000010.5.fits', SOURCEDIR + '\psc000010.6.fits']
RANGE = (0, 120)
m0collector = [None, None, None, [], [], [], []]

'''make header'''
HEADER =['ID']
#perImageHeaders = ['KronRad', 'separation', 'x', 'y', 'RA', 'DEC', 'KronMag', 'Angle', 'Ellipticity']
perImageHeaders = ['KronRad', 'separation', 'x', 'y','KronMag', 'Angle', 
                   'Ellipticity', 'RA', 'Host RA', 'DEC', 'Host Dec', 
                   'Discrepency', 'Z', 'Hecto Z']
for i in range(3,7):
    for val in perImageHeaders:
        HEADER.append(val + '_' + str(i))
                
                
# remove error file
if os.path.exists("errorfile.txt"):
  os.remove("errorfile.txt")
    
    
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
    plt.show()
    plt.close()

 
def extraction(filenames):
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
    global rad
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
    
    def errorProtocol(e, specified_file=None):
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
#            print "x"
        errorString = str(namecount) + '\n' + str(e) + " " + curFile + '\n'
        print errorString
        with open('errorfile.txt', 'a+') as errorfile:
            errorfile.write(errorString)
        if PLOT_ERR:
            plot(swappedData, objects, blacklist, bestCandidate, myVmin = 0, myVmax = 1000, 
                 myEventX = eventX, myEventY = eventY)
            plot(swappedData, objects, blacklist, bestCandidate, 
                 myEventX = eventX, myEventY = eventY)
        if ONLY_FLAG_ERRORS or e=='far' or e=='unlikely':
            return
        else:
            print "raising"
            raise
            
    ''' load event redshifts' dictionary '''
#TODO combine with plotRedshifts method
    zdict = {}
    zfile = open('new_ps1z.dat', 'r')
    zfile.readline() #get rid of header

    for line in zfile:
        parts = line.split() 
        eventId = parts[0][3:]
        redshift = float(parts[1])
        zdict[eventId] = redshift

    zfile.close()
    
    '''load event type dictionary'''
    typeDict = {}
    typefile = open('ps1confirmed_only_sne_without_outlier.txt', 'r')
    typefile.readline() #get rid of header
    for line in typefile:
        parts = line.split() 
        eventId = parts[0][3:]
        eventType = parts[1]
        typeDict[eventId] = eventType
        
    
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
    
#    if MAG_TEST_ALL or MAG_TEST_STDEV:
#        magdestfile = open(MAGFILE, "w+")
#        csvwriter = csv.writer(magdestfile)
    
    for filename in filenames:
        #clear bestCandidate to avoid bad plotting
        bestCandidate = None
        
        # to extract transient and filter number from filename of the form
        # */psc190369.6.fits
        dotSplit = filename.split('.')
        # transient identifier is everything after the last 'c' but before 
        # the second to last dot
        idNumString = dotSplit[-3].split('c')[-1] #used for getting hectospec data
        idNum  = int(idNumString)
        # filter number is between second to last dot and last dot
        filterNum = int(dotSplit[-2]) # filter of this specific image
        
        # if you are on a new event numbeer, 
        # fill the old row, write it, reset filter to 2, reset cache
        if lastIdNum != idNum:
            for imageDict in BAD_IMAGES:
                # these are for an event for which we found no host in center:
                errorProtocol('far', specified_file=str(lastIdNum) + '.' + \
                              str(imageDict['filterNum']) + ' no error info')
            GOOD_IMAGE = None
            BAD_IMAGES = []

            cache.extend([FILLER_VAL]*(1+4*len(perImageHeaders)-len(cache)))
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
        image_file = fits.open(filename)
        image_data = image_file[0].data
        
        #fix byte order
        swappedData = image_data.byteswap(inplace=True).newbyteorder()
        
        # subtracting out background
        bkg = sep.Background(swappedData)
        if SUBTRACT_BACKGROUND:
            swappedData = swappedData - bkg 
            
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
        
        '''Search SDSS for stars and blacklist them, new version'''
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
        
        image_center = SkyCoord(((maxRa + minRa) / 2.), ((maxDec + minDec) / 2.))
        circumscribing_diameter = SkyCoord(maxRa, maxDec).separation(SkyCoord(minRa, minDec))
        rad = (circumscribing_diameter/2.)
#TODO switch to rectangle rather than circle query 
        
        fields = ['ra','dec','type', 'mode', FILTERS[filterNum]]

        sdssTable = SDSS.query_region(coordinates=image_center, radius=rad, 
                                      photoobj_fields=fields)
        
        # objects to be ignored in host detection
        blacklist = set()
                           
        for i in range(len(objects)):
            # remove any 'nan' kronrads and blacklist the objects
            #check if ki is nan:
            if (not kronrad[i] == 0) and (not kronrad[i] < 0) and (not kronrad[i] > 0):
                # give it a valid kronrad, won't matter
                blacklist.add(i)
                
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

        # for collecting mags and fluxes to calculate the zero for this file
        colRealMags = []
        #colRealMags = Table([[]]*8, names=magNames)
        colFluxes = []
        for i in range(len(objects)):
            curRa = ra[i]
            curDec = dec[i]
            if sdssTable:    
                # for each object, iterate through table until you find primary object 
                #at that location
                for j in range(len(sdssTable)):
                    # it is a primary object not a duplicate
                    if sdssTable['mode'][j] == 1 \
                                and abs(sdssTable['ra'][j] - curRa)*u.deg < MINDIST \
                                and abs(sdssTable['dec'][j] - curDec)*u.deg < MINDIST:
                        # we have now found the primary sdss object matching object i's location                 
                        # if its a star, blacklist 
                        if sdssTable['type'][j] == 6: # is a star
                            blacklist.add(i)
                        
                        #if its a star brighter than 21, use its magnitude in zero point calculation
                        if sdssTable['type'][j] == 6 and sdssTable[FILTERS[filterNum]][j] < 21. \
                            and objects['x'][i] - 2.5 * kronrad[i] >= 0 \
                            and objects['x'][i] + 2.5 * kronrad[i] <= maxX \
                            and objects['y'][i] - 2.5 * kronrad[i] >= 0 \
                            and objects['y'][i] + 2.5 * kronrad[i] <= maxY:
                                colRealMags.append(sdssTable[FILTERS[filterNum]][j])
                                colFluxes.append(flux[i])
#TODO check that star is completely in frame before use                        
                        break   

        try:
            colFluxes = np.array(colFluxes)

#            magcache = [idNum, filterNum]
            m_0s = colRealMags + 2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS']))
            m_0s = m_0s[~np.isnan(m_0s)] #remove nans
            if m_0s.any(): 
                m_0 = np.median(m_0s)
                global m0collector
                m0collector[filterNum].append(m_0)

            #no SDSS stars to calculate zero point
            else: #or \

                    m_0 = FILTER_M0s[filterNum]['default']
                    colFluxes = np.array([])
                    colRealMags = np.array([])
            if m_0 > FILTER_M0s[filterNum]['upper'] or \
                m_0 < FILTER_M0s[filterNum]['lower']:  # m_0 is outlier
                    m_0 = FILTER_M0s[filterNum]['default']
            
            colMyMags = -2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS'])) + m_0
            magnitude = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + m_0
#TODO combine?
        
            
#            if MAG_TEST_STDEV:
#                csvwriter.writerow(magcache)
#                raise AssertionError
            if MAG_TEST_ALL:
                all_myMags[filterNum].extend(colMyMags[:])
                all_realMags[filterNum].extend(colRealMags[:])  
#                
#                for obj in range(len(colFluxes)):
#                    csvwriter.writerow([idNum, filterNum, colMyMags[obj], 
#                                        colRealMags[obj], 
#                                        colMyMags[obj] - colRealMags[obj]])    

            
                
#TODO remove
                for k in range(len(colMyMags)):
                    if abs(colMyMags[k] - colRealMags[k]) > 2:
                        print "outlier: " + filename
            
            '''Chance Coincidence Calculation'''
            #size is half light radius
            size, _flag = sep.flux_radius(swappedData, objects['x'], objects['y'], 
                              6.*objects['a'], 0.5, normflux=flux, subpix=5)
            # rough convert to arcsecs
# TODO refine?
            degPerPix = image_file[0].header['CDELT1']
            size = size * degPerPix *3600
              
            separation = objCoords.separation(eventCoords).arcsec

            # Observed number density of galaxies brighter than magnitude M (From Berger 2010)
            sigma = 10 ** (0.33 * (magnitude - 24) - 2.44) / (0.33 * np.log(10))
            # Effective radius
            R_effective = np.sqrt(np.abs(separation)**2 + 4 * np.abs(size) ** 2)
    
            # Probability of chance coincidence
            chanceCoincidence = 1 - np.exp(-np.pi * R_effective ** 2 * sigma)
            # exclude any blacklisted
            for i in blacklist:
                chanceCoincidence[i] = 1
            bestCandidate = np.argmin(chanceCoincidence)
                
            ellipticity = 1 - (objects['b'][bestCandidate]/objects['a'][bestCandidate])
    
            image_file.close()
#TODO remove repitition of calculations between filter nums    
        
            '''Get "real" host location and redshifts'''
            with open('hosts.dat', 'r') as hostsfile:
                hostsData = hostsfile.read()
            # convert to dictionary from dictionary string
            hostsData = ast.literal_eval(hostsData)
            hostRa = hostsData[idNumString]['host_ra']
            hostDec = hostsData[idNumString]['host_dec']
            try:
                hostCoords = SkyCoord(hostRa, hostDec, unit=(u.hourangle, u.deg))
                #convert to decimal deg:
                hostRa = hostCoords.ra.deg
                hostDec = hostCoords.dec.deg            
                offby = hostCoords.separation(objCoords[bestCandidate]).arcsec\
                        if hostRa else None
                hectoZ = hostsData[idNumString]['redshift']
            except ValueError: 
                if len(hostRa) > 12: # hecto gave multiple host galaxies
                    hostRa = "multiple"
                    hostDec = "multiple"
                    offby = None
                    hectoZ = "multiple"
                else:
                    raise
            
            if PLOT_ALL:
                plot(swappedData, objects, blacklist, bestCandidate,
                     myEventX = eventX, myEventY = eventY)
                plot(swappedData, objects, blacklist, bestCandidate, 
                     myVmin = 0, myVmax = 1000, myEventX = eventX, myEventY = eventY)
            try:
                z = zdict[idNumString]
            except KeyError: # no redshift for this event:
                z = None
              
            finalProperties = [kronrad[bestCandidate], separation[bestCandidate], 
                    objects['x'][bestCandidate] - eventX, 
                    objects['y'][bestCandidate] - eventY,
                    magnitude[bestCandidate], objects['theta'][bestCandidate], 
                    ellipticity, ra[bestCandidate], hostRa, dec[bestCandidate], 
                    hostDec, offby, z, hectoZ]

            if separation[bestCandidate] > CHECK_DISTANCE:
                #TODO remove any unnecessary
                badImageDict = {'filterNum': filterNum, 'objects': objects, 
                                'kronrad': kronrad, 'separation': separation,
                                'magnitude': magnitude, 'objCoords': objCoords,
                                'm_0': m_0}
                BAD_IMAGES.append(badImageDict)
            else:
                newFluxParams =  (swappedData, 
                               objects['x'][bestCandidate],
                               objects['y'][bestCandidate],
                               objects['a'][bestCandidate],
                               objects['b'][bestCandidate], 
                               objects['theta'][bestCandidate], 2.5*kronrad[bestCandidate])
                
                GOOD_IMAGE = (finalProperties, newFluxParams)
                
            cache.extend(finalProperties)            
            
            # If possible, repair cache from bad images, whether this or previous filters
            if GOOD_IMAGE and BAD_IMAGES: 
                print "correcting"
                finalProperties, newFluxParams = GOOD_IMAGE
                #TODO check indices are correct
                goodCoords = SkyCoord(finalProperties[8], 
                                      finalProperties[10], unit=u.deg) #(ra,dec)
                for imageDict in BAD_IMAGES:
                    fixed = False
                    for k in range(len(badImageDict['objCoords'])):
                        # look for image 
                        if badImageDict['objCoords'][k].separation(goodCoords) < MINDIST:
                            newBestCandidate = k
                            newEllipticity = 1 - (imageDict['objects']['b'][newBestCandidate]/
                                                imageDict['objects']['a'][newBestCandidate])
                            newOffby = hostCoords.separation(imageDict['objCoords'][newBestCandidate])
                            newProperties = [imageDict['kronrad'][newBestCandidate], 
                                            imageDict['separation'][newBestCandidate],
                                            imageDict['objects']['x'][newBestCandidate] - eventX,
                                            imageDict['objects']['y'][newBestCandidate] - eventY,
                                            imageDict['magnitude'][newBestCandidate],
                                            imageDict['objects']['theta'][newBestCandidate],
                                            newEllipticity,
                                            imageDict['objCoords'][newBestCandidate].ra.deg,
                                            hostRa, 
                                            imageDict['objCoords'][newBestCandidate].dec.deg,
                                            hostDec, newOffby, z, hectoZ]
                            fixed = True
                            break
                    if not fixed: # no objects in center detected in bad image, use data from good image
                        # but update magnitude
                        newFlux, _fluxerr, _flag = sep.sum_ellipse(*newFluxParams, subpix=1)
                        newMagnitude = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + imageDict['m_0']
                        newProperties = finalProperties[:] #copy list
                        newProperties[5] = newMagnitude
                    #TODO fix out of bounds possibility
                    numColumnsOnLeft = 1 + (imageDict['filterNum'] - 3)*len(perImageHeaders)
                    cache = cache[:numColumnsOnLeft] + newProperties\
                            + cache[numColumnsOnLeft + len(perImageHeaders):]
                    BAD_IMAGES.remove(imageDict)    
            
            if chanceCoincidence[bestCandidate] > LIKELIHOOD_THRESH:
                errorProtocol("unlikely")
                    
            ''' end getData block'''    
            #collect data for redshift plot
            if PLOT_REDSHIFTS and filterNum == 6:    
                thisType = typeDict[idNumString]
                all_redshifts[thisType].append(cache[-2])
                all_kronMags[thisType].append(cache[-10])
                all_kronRads[thisType].append(cache[-14])
                    
            
            lastIdNum = idNum
            lastFilter = filterNum
             
#TODO remove
        except AssertionError:
            print "here"
            if MAG_TEST_ALL or MAG_TEST_STDEV:
                continue
            else:
                raise
        except Exception as e:
            errorProtocol(e)
    if PLOT_REDSHIFTS:
        #TYPE_COLORS = {'SNIIn':'co', 'SNIa':'ro', 'SNII': 'bo', 'SNIbc':'go', 'SLSNe': 'mo'}
        
        # plot magnitude vs. redshift
        for snType in TYPES:
            plt.plot(all_redshifts[snType], all_kronMags[snType], TYPE_COLORS[snType])
          
        plt.xlabel('Redshifts')
        plt.ylabel('Magnitudes')
        plt.savefig('redshifts_vs_magnitudes', dpi=150)
        plt.show()
        plt.close()
        
        areas = {}
        # plot magnitude vs. area
        for snType in TYPES:
            all_kronMags[snType] = np.array(all_kronMags[snType])
            all_kronRads[snType] = np.array(all_kronRads[snType])
            areas[snType] = np.pi * all_kronRads[snType]**2
            plt.plot(all_redshifts[snType], areas[snType], TYPE_COLORS[snType])
          
        plt.xlabel('Redshifts')
        plt.ylabel('Surface Areas')
        plt.savefig('redshifts_vs_areas', dpi=150)
        plt.show()
        plt.close()
        
        surface_brightness = {}
        # plot magnitude vs. surface brightness
        for snType in TYPES:
            surface_brightness[snType] = all_kronMags[snType]/areas[snType]
            plt.plot(all_redshifts[snType], surface_brightness[snType], TYPE_COLORS[snType])
          
        plt.xlabel('Redshifts')
        plt.ylabel('Surface Brightness')
        plt.savefig('redshifts_vs_brightness', dpi=150)
        plt.show()
        plt.close()
        
    if MAG_TEST_ALL or MAG_TEST_STDEV:
#        all_myMags2 = np.array(all_myMags)
#        all_realMags2 = np.array(all_realMags)
#        #magdestfile.close()    
#        plt.plot(all_myMags2, all_realMags2 - all_myMags2, 'bo')
#        plt.xlabel('My magnitudes')
#        plt.ylabel('Difference')#'SDSS magnitudes')
#        plt.savefig('magplot_mymags_vs_dif_leaving_outliers_newDefault', dpi=150)
#        plt.show()
#        plt.close()
#        
#        plt.plot(all_realMags2, all_realMags2 - all_myMags2, 'bo')
#        plt.xlabel('SDSS magnitudes')
#        plt.ylabel('Difference')#'SDSS magnitudes')
#        plt.savefig('magplot_sdssmags_vs_dif_leaving_outliers_newDefault', dpi=150)
#        plt.show()
#        plt.close()
#        
#        plt.plot(all_myMags2, all_realMags2, 'bo')
#        plt.xlabel('My magnitudes')
#        plt.ylabel('SDSS Magnitudes')#'SDSS magnitudes')
#        plt.savefig('magplot_sdssmags_vs_mymags_leaving_outliers_newDefault', dpi=150)
#        plt.show()
#        plt.close()
#        
#        # plot again w/o median points
#        all_myMags3 = []
#        all_realMags3 = []
#        for l in range(len(all_myMags2)):
#            if all_myMags2[l] != all_realMags2[l]:
#                all_myMags3.append(all_myMags2[l])
#                all_realMags3.append(all_realMags2[l])
#       
#        all_myMags4 = np.array(all_myMags3)
#        all_realMags4 = np.array(all_realMags3)
        global all_realMagsFiltered
        global all_myMagsFiltered
        all_realMagsFiltered = [[], [], [], [], [], [], []]
        all_myMagsFiltered = [[], [], [], [], [], [], []]
        for x in range(3,7):
            for l in range(len(all_myMags[x])):
                if all_myMags[x][l] != all_realMags[x][l]:
                    all_myMagsFiltered[x].append(all_myMags[x][l])
                    all_realMagsFiltered[x].append(all_realMags[x][l])
                    
            all_realMagsFiltered[x] = np.array(all_realMagsFiltered[x])
            all_myMagsFiltered[x] = np.array(all_myMagsFiltered[x])
            print x
            
        plt.plot(all_myMagsFiltered[3], all_realMagsFiltered[3] - all_myMagsFiltered[3], 'bo',
                 all_myMagsFiltered[4], all_realMagsFiltered[4] - all_myMagsFiltered[4], 'go',
                 all_myMagsFiltered[5], all_realMagsFiltered[5] - all_myMagsFiltered[5], 'ro',
                 all_myMagsFiltered[6], all_realMagsFiltered[6] - all_myMagsFiltered[6], 'yo',)
        plt.xlabel('My magnitudes')
        plt.ylabel('Difference')#'SDSS magnitudes')
        plt.savefig('deletable', dpi=150)
        plt.show()
        plt.close()
            #magplot_mymags_vs_dif_nonzero_outliersCoerced_noGalaxies_noDim
#            plt.plot(all_realMagsFiltered, all_realMagsFiltered - all_myMagsFiltered, 'bo')
#            plt.xlabel('SDSS magnitudes')
#            plt.ylabel('Difference')#'SDSS magnitudes')
#            plt.savefig('magplot_sdssmags_vs_dif_nonzero_leaving_outliers_newDefault'+str(x), dpi=150)
#            plt.show()
#            plt.close()
#            
#            plt.plot(all_myMagsFiltered, all_realMagsFiltered, 'bo')
#            plt.xlabel('My magnitudes')
#            plt.ylabel('SDSS Magnitudes')#'SDSS magnitudes')
#            plt.savefig('magplot_sdssmags_vs_mymags_nonzero_leaving_outliers_newDefault'+str(x), dpi=150)
#            plt.show()
#            plt.close()
#        
    # write last data
    cache.extend([FILLER_VAL]*(1+4*len(perImageHeaders)-len(cache)))
    if WRITE_CSV:
        csvwriter.writerow(cache)
        destfile.close()
    if PRINT_DATA:
        print cache
   
def main():
#TODO remove
    import time
    start = time.time()
    #figure out which files to use based on value specified at top
    if FILES == 'all' or FILES =='range' or FILES =='new random':
        filenames = sorted(glob.glob(SOURCEDIR + '/psc*.[3-6].fits'))
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
    print end - start

if __name__ == "__main__":
     main()
