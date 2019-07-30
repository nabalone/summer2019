# -*- coding: utf-8 -*-

#For debugging:
PLOT = False
CIRCLE_OBJS = False #only set to true if PLOT
CIRCLE_HOST = False



"""
https://docs.python.org/3/library/csv.html
# http://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html
# http://docs.astropy.org/en/stable/io/fits/


find todos
remove background removal?
"""
import numpy as np
#from matplotlib.patches import Ellipse
import os
import glob
import csv
import sep
#import matplotlib.pyplot as plt
#from astropy.visualization import astropy_mpl_style
#plt.style.use(astropy_mpl_style)
#from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd() + "\galaxiestest"
FILLER_VAL = None
THRESHOLD = 3
MINAREA = 5
MAG_ZERO =25 # for calculating magnitude from flux
DEBLEND_CONT = 0.05 # for sep.extract. 1.0 to turn of deblending (works best),
 # 0.005 is default
SUBTRACT_BACKGROUND = True

#TODO remove, for debugging
errors = []

# make header
HEADER =['ID']
perImageHeaders = ['KronRad', 'separation', 'x', 'y', 'RA', 'DEC', 'KronMag', 'Angle', 'Ellipticity']
for i in range(3,7):
    for val in perImageHeaders:
        HEADER.append(val + '_' + str(i))


filenames = sorted(glob.glob(SOURCEDIR + '/psc*'))

#helper function to extract the data from a file
def getData(filename, threshold=THRESHOLD):
#TODO does this need to be recursive? if so move the following out, if not maybe
    #get rid of helper fxn
    image_file = fits.open(filename)
    image_data = image_file[0].data

    # https://www.aanda.org/index.php?option=com_article&access=doi&doi=
    # 10.1051/0004-6361/201015362&Itemid=129#FD14
    
    #fix byte order
    swappedData = image_data.byteswap(inplace=True).newbyteorder()
    
    # subtracting out background
    bkg = sep.Background(swappedData)
    if SUBTRACT_BACKGROUND:
        swappedData = swappedData - bkg 
    
    objects = sep.extract(swappedData, threshold, err=bkg.globalrms, 
                          minarea = MINAREA, deblend_cont = DEBLEND_CONT)
    
 #TODO if not objects:
    
    
    '''Find most likely host galaxy'''
    # how to calculate kron radius and flux from 
    # https://sep.readthedocs.io/en/v1.0.x/apertures.html
    kronrad, _krflag = sep.kron_radius(swappedData, objects['x'], objects['y'], 
                                       objects['a'], objects['b'], 
                                       objects['theta'], 6.0)
    print len(objects)
    # list of bad objects to be ignored
    blacklist = []
    # remove any nan kronrads and blacklist the objects
    for i in range(len(objects)):
        #check if ki is nan:
        if (not kronrad[i] == 0) and (not kronrad[i] < 0) and (not kronrad[i] > 0):
            # remove that object's column
            kronrad[i] = 0.1
            blacklist.append(i)
    flux, _fluxerr, _flag = sep.sum_ellipse(swappedData, objects['x'], 
                                            objects['y'], objects['a'], 
                                            objects['b'], objects['theta'], 
                                            2.5*kronrad, subpix=1)
    magnitude = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + MAG_ZERO 
    
    #size is half light radius
    size, _flag = sep.flux_radius(swappedData, objects['x'], objects['y'], 
                      6.*objects['a'], 0.5, normflux=flux, subpix=5)
        
    separation = np.sqrt(((objects['x']-120)**2) + ((objects['y']-120)**2))
    
    # Observed number density of galaxies brighter than magnitude M (From Berger 2010)
    sigma = 10 ** (0.33 * (magnitude - 24) - 2.44) / (0.33 * np.log(10))

    # Effective radius
    R_effective = np.sqrt(np.abs(separation) + 4 * np.abs(size) ** 2)

    # Probability of chance coincidence
    chanceCoincidence = 1 - np.exp(-np.pi * R_effective ** 2 * sigma)
#TODO check reasonable chance coincidence 
    # exclude any blacklisted
    for i in blacklist:
        chanceCoincidence[i] = 1
    bestCandidate = np.argmin(chanceCoincidence)
    
    if separation[bestCandidate] > 20:
        #print("far: " + filename + " "+ str(separation[bestCandidate]))
        pass
    ellipticity = 1- (objects['b'][bestCandidate]/objects['a'][bestCandidate])
    
    '''  
    -2.5 * np.log10(flux) + MAG_ZERO
    f.write(filename)    
    f.write('\n'+ str(magnitude[bestCandidate]))
    f.write('\n')
    mag2 = -2.5 * np.log10(flux/float(image_file[0].header['MJD-OBS'])) + MAG_ZERO
    f.write(str(mag2[bestCandidate]))
    f.write('\n\n\n\n')
    
'''
#TODO make sure directions are correct
    centerRa = image_file[0].header['CRVAL1']
    centerDec = image_file[0].header['CRVAL2']
    ra = centerRa-image_file[0].header['CDELT1'] * objects['x'][bestCandidate]
    dec = centerDec-image_file[0].header['CDELT2']*objects['y'][bestCandidate]
    image_file.close()
    return [kronrad[bestCandidate], separation[bestCandidate], 
            objects['x'][bestCandidate], objects['y'][bestCandidate],
            ra, dec, magnitude[bestCandidate],
            objects['theta'][bestCandidate], ellipticity]

def main():
    with open(DESTDIR + "\galaxiesdata3.csv", "w+") as destfile:
        csvwriter = csv.writer(destfile)
        csvwriter.writerow(HEADER)
        
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
                csvwriter.writerow(cache)
                cache = [idNum]
                lastFilter = 2
                
            # fill any missing data in current row
            for i in range(filterNum - lastFilter - 1):
                cache.extend([FILLER_VAL]*len(perImageHeaders))
#TODO remove
            try:
                cache.extend(getData(filename))
            except:
                #pass
                print("error: " + filename)
                #raise
                errors.append(filename)
        
            lastIdNum = idNum
            lastFilter = filterNum
        
        # write last data
        cache.extend([FILLER_VAL]*(29-len(cache)))
        csvwriter.writerow(cache)
        cache = [idNum]
        lastFilter = 2

        
if __name__ == "__main__":
    #main()
    print(getData('C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc470240.3.fits'))
    #print(getData('C:/Users/Faith/Desktop/noey2019summer/ps1hosts\psc000010.6.fits'))
'''  
f = open('magstest2.txt', 'w+')
for i in range(80,100):
    print filenames[i]
    getData(filenames[i], f)  
    print("done")
f.close()

            # check if this is the last filter(6) of transient, if so write the 
            #cache to create line for this trasient. reset cache 
            if filterNum == '6':
                if len(cache) != 29:
                    raise Exception("Incorrect row length")
                csvwriter.writerow(cache)
                cache = []
                
 '''   
