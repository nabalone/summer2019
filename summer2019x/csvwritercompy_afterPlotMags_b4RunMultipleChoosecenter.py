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
import shutil
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

SINGLE_QUERY = True  
SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd() + "\galaxiestest2"
FILLER_VAL = None
THRESHOLD = 1.5
MINAREA = 5
#TODO remove?
MAG_ZERO =25 # for calculating magnitude from flux
DEBLEND_CONT = 0.05 # for sep.extract. 1.0 to turn of deblending (works best),
 # 0.005 is default
SUBTRACT_BACKGROUND = True
MINDIST = 0.0005*u.deg #dist. an sdss object must be within to identify as
FILTERS = [None, None, None, 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z']
LIKELIHOOD_THRESH = 0.9

#For figuring calibrating magnitudes
#magNames = ['modelMag_z', 'modelMag_i', 'modelMag_r', 'modelMag_g',
#            'psfMag_z', 'psfMag_i', 'psfMag_r', 'psfMag_g']

#For debugging:
PLOT_REDSHIFTS_ONLY = False
WRITE_CSV = None# "\galaxiesdata3.csv" # filename to write to or None
MAG_TEST_ALL = True
MAG_TEST_STDEV = False
MAGFILE = 'newmagtest.csv'

#Do not write csv an mag test simultaneously because there is only 1 writer
#TODO separate the csvwriters
if MAG_TEST_ALL or MAG_TEST_STDEV:
    WRITE_CSV = False
PRINT_DATA = False
CHECK_DISTANCE = 20 #print all files with most likely host farther than this int
PLOT_ALL = False
PLOT_ERR =  False #plots only files that give errors or low probability
PLOT_DIR = os.getcwd() + '/plotall' # where to put plot images
# if plot_dir already exists, clear it
#shutil.rmtree(PLOT_DIR)
ONLY_FLAG_ERRORS = False # catch errors, print filename, move on
FILES = 'range' #options are 'all', 'preset random', 'new random', 'range', 'specified', 'nonsquare
SPECIFIED = [SOURCEDIR + '\psc000010.3.fits']
#[SOURCEDIR + '\psc000010.3.fits', SOURCEDIR + '\psc000010.4.fits', SOURCEDIR + '\psc000010.5.fits', SOURCEDIR + '\psc000010.6.fits']
RANGE = (0, 60)


# make header
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

 #TODO remove
''' My stupid workaround for debugging, so I can access the variables from the console '''
objects=None
bestCandidate=None
magnitudes=None
sdssTable=None
chanceCoincidence=None
separation=None
filename=None
blacklist=None
zs = [] #redshifts
ms = [] #magnitudes with indices corresponding to zs
image_file=None
swappedData = None
rad = None
zdict = {}
hostsData = None
all_realMags = []
all_myMags = []
    
    
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
    

def plotRedshifts():
    ''' My stupid workaround for debugging, so I can access the variables from the console'''
    global zs #redshifts
    global ms #magnitudes with indices corresponding to zs
    zfile = open('new_ps1z.dat', 'r')
    zfile.readline() #get rid of header
    # iterate through list of redshifts, find corresponding .fits file and
    # extract magnitude
    for i in range(50):
        line = zfile.readline()
        parts = line[:-1].split('\t') #strip line ending and split by tabs
        eventId = parts[0][3:]
        redshift = float(parts[1])
        filename = SOURCEDIR + '/psc' + eventId + '.6.fits'
        try:
            mag = extraction([filename], mag_only=True)
            ms.append(mag)
            zs.append(redshift)
        except IOError as e: #probably the .fits file was missing
            print(e)
    zfile.close()
    
    plt.plot(zs, ms, 'bo')
    plt.savefig('redshift_vs_magnitude', dpi=150)
 
def extraction(filenames, mag_only=False):
    ''' My stupid workaround for debugging, so I can access the variables from the console'''
    global objects
    global bestCandidate
    global magnitudes
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
    global zs #redshifts
    global ms #magnitudes with indices corresponding to zs
    global zdict
    
    def errorProtocol(e):
        print("error: " + filename + ' ' + str(e))
        with open('errorfile.txt', 'a+') as errorfile:
            errorfile.write(str(namecount))
            errorfile.write("error: " + filename + ' ' + str(e) + '\n')
        if PLOT_ERR:
            plot(swappedData, objects, blacklist, bestCandidate, myVmin = 0, myVmax = 1000, 
                 myEventX = eventX, myEventY = eventY)
            plot(swappedData, objects, blacklist, bestCandidate, 
                 myEventX = eventX, myEventY = eventY)
        if ONLY_FLAG_ERRORS:
            return
        else:
            print "raising"
            raise
            
    ''' load event redshifts' dictionary '''
#TODO combine with plotRedshifts method

    zdict = {}
    zfile = open('new_ps1z.dat', 'r')
    zfile.readline() #get rid of header
    # iterate through list of redshifts, find corresponding .fits file and
    # extract magnitude
    #for i in range(50):
        #line = zfile.readline()
    for line in zfile:
        parts = line[:-1].split('\t') #strip line ending and split by tabs
        eventId = parts[0][3:]
        redshift = float(parts[1])
        zdict[eventId] = redshift
#                filename = SOURCEDIR + '/psc' + eventId + '.6.fits'
#                try:
#                    mag = extraction([filename], mag_only=True)
#                    ms.append(mag)
#                    zs.append(redshift)
#                except IOError as e: #probably the .fits file was missing
#                    print(e)
    zfile.close()

#    plt.plot(zs, ms, 'bo')
#    plt.savefig('redshift_vs_magnitude', dpi=150)
        
    
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
    
    if MAG_TEST_ALL or MAG_TEST_STDEV:
        magdestfile = open(MAGFILE, "w+")
        csvwriter = csv.writer(magdestfile)
    
    for filename in filenames:


        # to extract transient and filter number from filename of the form
        # */psc190369.6.fits
        dotSplit = filename.split('.')
        # transient identifier is everything after the last 'c' but before 
        # the second to last dot
        idNumString = dotSplit[-3].split('c')[-1] #used for getting hectospec data
#TODO is this necessary?
        idNum  = int(idNumString)
        # filter number is between second to last dot and last dot
        filterNum = int(dotSplit[-2]) # filter of this specific image
        
        # if you are on a new event numbeer, 
        # fill the old row, write it, reset filter to 2, reset cache
        if lastIdNum != idNum:
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
        if MAG_TEST_ALL or MAG_TEST_STDEV:
            # for collecting mags and fluxes to calculate the zero for this file
            colRealMags = []
            #colRealMags = Table([[]]*8, names=magNames)
            colFluxes = []
        
        
        '''Search SDSS for stars and blacklist them, new version'''
# TODO remove SINGLE_QUERY VARIABLE
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
        #WHY WON'T IT LET ME SPECIFY WIDTH???
        # use circle circumscribing image until I can figure out how to use width
        #rad = (np.sqrt(2) * imageWidth / 2.) * u.deg
        fields = ['ra','dec','type', 'mode'] #filter
        if MAG_TEST_ALL or MAG_TEST_STDEV:
            fields.append(FILTERS[filterNum])

        sdssTable = SDSS.query_region(coordinates=image_center, radius=rad, 
                                      photoobj_fields=fields)
                                         #width=imageWidth*u.deg,height=imageHeight)                                   
        for i in range(len(objects)):
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
        except Exception as e:
            errorProtocol(e)
                
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
                        # collect its data then break
                        
                        if (MAG_TEST_ALL or MAG_TEST_STDEV) and sdssTable['type'][j] == 6:
                            colRealMags.append(sdssTable[FILTERS[filterNum]][j])
#                            curMags = []
#                            for n in magNames:
#                                curMags.append(sdssTable[n][j])
#                            # flux[x] corresponds to the xth row of colRealMags
#                            colRealMags.add_row(curMags)
                            colFluxes.append(flux[i])
                       
                        # if its a star, blacklist 
                        if sdssTable['type'][j] == 6: # is a star
                            blacklist.add(i)
                            
                        break

        try:
            if MAG_TEST_ALL or MAG_TEST_STDEV:
                global colFluxes
                colFluxes = np.array(colFluxes)
                
                #collect m_0s
                #final_m_0[x] corresponds to the filter magNames[x]
                #final_m_0s = {}
                #colMyMags = {}
                global m_0s
                global m_0
                magcache = [idNum, filterNum]
                m_0s = colRealMags + 2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS']))
                m_0s = m_0s[~np.isnan(m_0s)] #remove nans
                m_0 = np.median(m_0s)
            
                colMyMags = -2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS'])) + m_0
                
                # alternate calculation, I have verified it matches
#                alt_offsetMags = - 2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS']))
#                alt_difs = colRealMags - alt_offsetMags
#                alt_medDif = np.median(alt_difs)
#                alt_colMyMags = alt_offsetMags + alt_medDif
                
                
#               for filt in magNames:
                    #format filter column into an ndarray
    #TODO clean up line below
#                    colFilterMags = np.array(str(colRealMags[filt]).split()[2:], dtype = 'float64')
#                    m_0s = colFilterMags + 2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS']))
#    
#                    magcache.extend([filt, np.std(m_0s)])
#                    
#                    m_0 = np.median(m_0s)
#                    final_m_0s[filt] = m_0
#                    colMyMags = -2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS'])) + m_0
                #TODO remove!! for magstd testing!
            if MAG_TEST_STDEV:
                csvwriter.writerow(magcache)
                raise AssertionError
            elif MAG_TEST_ALL:
                global all_mymags
                global all_realmags
                all_myMags.extend(colMyMags)
                all_realMags.extend(colRealMags)
                for obj in range(len(colFluxes)):
#                    magCache = [idNum, filterNum]
#                    for filter in magNames:
#                        magCache.extend([None, filter[-1], final_m_0s[filter], colMyMags[filter][obj], 
#                                        colRealMags[filter][obj], colMyMags[filter][obj] - colRealMags[filter][obj]])
                    csvwriter.writerow([idNum, filterNum, colMyMags[obj], 
                                        colRealMags[obj], colMyMags[obj] - colRealMags[obj]])    
                #csvwriter.writerow(magCache)
                #csvwriter.writerow([None]*50)
            # colRealMags = np.array(colRealMags)
            # m_0s = colRealMags + 2.5 * np.log10(colFluxes)
            # m_0 = np.median(m_0s)
#TODO change MAG_ZERO to an m_0
            global magnitude
            magnitude = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + MAG_ZERO
            if MAG_TEST_ALL or MAG_TEST_STDEV:
                colMyMags = -2.5 * np.log10(colFluxes/float(image_file[0].header['MJD-OBS'])) + m_0
            global size
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
            global R_effective
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
                if separation[bestCandidate] > CHECK_DISTANCE \
                    or chanceCoincidence[bestCandidate] > LIKELIHOOD_THRESH:
                    print("far: " + filename + " "+ str(separation[bestCandidate]))
                    with open('errorfile.txt', 'a+') as errorfile:
                        errorfile.write(filename + 
                                        " dist: " + str(separation[bestCandidate]) + \
                                        " chanceCoincidence: " + str(chanceCoincidence[bestCandidate])
                                        + '\n')
                    if PLOT_ERR:
                        plot(swappedData, objects, blacklist, bestCandidate, myVmin = 0, myVmax = 1000, 
                             myEventX = eventX, myEventY = eventY)
                        plot(swappedData, objects, blacklist, bestCandidate, 
                             myEventX = eventX, myEventY = eventY)
            ellipticity = 1 - (objects['b'][bestCandidate]/objects['a'][bestCandidate])
    
            image_file.close()
            
            '''Get "real" host location and redshifts'''
            global hostsData
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
                offby = hostCoords.separation(eventCoords).deg if hostRa else None
                hectoZ = hostsData[idNumString]['redshift']
            except ValueError: 
                if len(hostRa) > 12: # hecto gave multiple host galaxies
                    hostRa = "multiple"
                    hostDec = "multiple"
                    hectoZ = "multiple"
                else:
                    raise

            
                
            cache.extend([kronrad[bestCandidate], separation[bestCandidate], 
                    objects['x'][bestCandidate] - eventX, objects['y'][bestCandidate] - eventY,
                    magnitude[bestCandidate], objects['theta'][bestCandidate], 
                    ellipticity, ra[bestCandidate], hostRa, dec[bestCandidate], 
                    hostDec, offby, zdict[idNumString], hectoZ])
                    
            ''' end getData block'''            
            
            lastIdNum = idNum
            lastFilter = filterNum
            
    #        print filename
    #        print objCoords.to_string('hmsdms', sep=':')
    #        print magnitude
    #        print ''
    
            if mag_only:
                return magnitude[bestCandidate]
            
            if PLOT_ALL:
                plot(swappedData, objects, blacklist, bestCandidate,
                     myEventX = eventX, myEventY = eventY)
                plot(swappedData, objects, blacklist, bestCandidate, 
                     myVmin = 0, myVmax = 1000, myEventX = eventX, myEventY = eventY)
             
            #print chanceCoincidence
#TODO remove
        except AssertionError:
            print "here"
            if MAG_TEST_ALL or MAG_TEST_STDEV:
                continue
            else:
                raise
        except Exception as e:
            errorProtocol(e)
    if MAG_TEST_ALL or MAG_TEST_STDEV:
        all_myMags2 = np.array(all_myMags)
        all_realMags2 = np.array(all_realMags)
        magdestfile.close()    
        plt.plot(all_myMags2, all_realMags2 - all_myMags2, 'bo')
        plt.xlabel('My magnitudes')
        plt.ylabel('Difference')#'SDSS magnitudes')
        plt.show()
        plt.savefig('magplot_mymags_vs_dif', dpi=150)
        
        plt.plot(all_realMags2, all_realMags2 - all_myMags2, 'bo')
        plt.xlabel('SDSS magnitudes')
        plt.ylabel('Difference')#'SDSS magnitudes')
        plt.show()
        plt.savefig('magplot_sdssmags_vs_dif', dpi=150)
        
        plt.plot(all_myMags2, all_realMags2, 'bo')
        plt.xlabel('My magnitudes')
        plt.ylabel('SDSS Magnitudes')#'SDSS magnitudes')
        plt.show()
        plt.savefig('magplot_sdssmags_vs_dif', dpi=150)
        
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
    if PLOT_REDSHIFTS_ONLY:
        plotRedshifts()
    #figure out which files to use based on value specified at top
    elif FILES == 'all' or FILES =='range' or FILES =='new random':
        filenames = sorted(glob.glob(SOURCEDIR + '/psc*.[3-6].fits'))
        if FILES == 'range':
            extraction(filenames[RANGE[0]:RANGE[1]])
        elif FILES == 'new random':
            randFiles = []
            for i in range(10):
                randFiles.append(filenames[int(len(filenames)*random.random())])
            extraction(randFiles)
        elif FILES =='all':
            extraction(filenames)
    elif FILES == 'preset random':
        from presetrandomfiles import fileset
        extraction(fileset[:25])
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

'''#TODO update with mag checking or delete    
        else:
            for i in range(len(objects)):
                sdssTable = SDSS.query_region(coordinates=objCoords[i], radius=MINDIST, 
                                      photoobj_fields=['ra','dec','type', 'mode'])
        
                #if table of objects at that location exists and 
                # a star and no galaxies are there. star is type 6, galaxy is type 3
                #if sdssTable and 6 in sdssTable['type'] and not 3 in sdssTable['type']:
                # nevermind, actually if the primary object (mode=1) is star then blacklist
                
                if sdssTable:    
                    for j in range(len(sdssTable)):
                        if sdssTable['mode'][j] == 1 and sdssTable['type'][j] == 6:
                            blacklist.add(i)
                            break
                # remove any 'nan' kronrads and blacklist the objects
                #check if ki is nan:
                if (not kronrad[i] == 0) and (not kronrad[i] < 0) and (not kronrad[i] > 0):
                    # remove that object's column
                    kronrad[i] = 0.1
                    blacklist.add(i)
                    
                    '''
