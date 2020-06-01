# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:02:48 2019

@author: Faith
"""

# Note: filters 3,4,5,6 are g,r,i,z respectively
# Types 0,1,2,3,4 are types Ia, Ibc, II, IIn, and superluminous respectively

import numpy as np
import pandas as pd
import os
import glob
import sys
#import matplotlib.pyplot as plt
#from matplotlib.patches import RegularPolygon
#from astropy.visualization import astropy_mpl_style
#plt.style.use(astropy_mpl_style)
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import sep

PROJ_HOME = os.environ['DATA_SRCDIR']
sys.path.append(PROJ_HOME)
#from src.random_forest_classifier import chooseAll

#supress fits warnings
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


SOURCEDIR = os.getcwd() 
BASEDIR = PROJ_HOME + "/src/"
OUTPUT_DIR = BASEDIR + "outputs/" #where 'all_masks.npz' is located
PIXDIR = BASEDIR + "all_fits/"
CSVFILE = BASEDIR + 'outputs/galaxiesdata.csv' 

filenames = sorted(glob.glob(PIXDIR + 'psc*.[3].fits'))
x_train = []
y_train = []
x_test = []
y_test = []

'''load event type dictionary'''
typeDict = {}
typefile = open(BASEDIR + 'ps1confirmed_added.txt', 'r')
typefile.readline() #get rid of header
for line in typefile:
    parts = line.split()
    eventId = parts[0][3:]
    eventType = parts[1]
    typeDict[eventId] = eventType
typefile.close()

zdict = {}
zfile = open(BASEDIR + 'new_ps1z.dat', 'r')
zfile.readline() #get rid of header

for line in zfile:
    parts = line.split()
    eventId = parts[0][3:]
    redshift = float(parts[1])
    zdict[eventId] = redshift
zfile.close()

intDict= {'SNIa':0, 'SNIbc':1, 'SNII':2, 'SNIIn':3, 'SLSNe':4}

'''load event location dictionary'''
db = pd.read_table(BASEDIR + 'alertstable_v3',sep=None,index_col = False,
               engine='python')
db2 = pd.read_table(BASEDIR + 'alertstable_v3.lasthalf',sep=None,index_col = False,
                engine='python')
db = db.append(db2,ignore_index=True)

def size(axis, image_data):
    return int(np.shape(image_data)[axis])

def pad(axis, side, image_data):
    width = 240 - size(axis, image_data)
    if axis == 0:
        zeros = np.zeros((width, size(1, image_data)))
        if side == 'below':
            return np.vstack((zeros, image_data))
        elif side == 'above':
            return np.vstack((image_data, zeros))
    else:
        zeros = np.zeros((size(0, image_data), width))
        if side == 'below':
            return np.hstack((zeros, image_data))
        elif side == 'above':
            return np.hstack((image_data, zeros))
        

#def plot(data, myEventX, myEventY, myVmin=None, myVmax=None):
#    fig, ax = plt.subplots()
#
#    if myVmax == None or myVmin == None:
#        _im = ax.imshow(data, interpolation='nearest', cmap='gray')
#    else:
#        _im = ax.imshow(data, interpolation='nearest', cmap='gray',
#                        vmin = myVmin, vmax = myVmax)
#        
#    # triangle on event location
#    p = RegularPolygon((int(myEventX), int(myEventY)), 3, radius=3)
#    p.set_edgecolor('purple')
#    p.set_facecolor('purple')
#    ax.add_artist(p)
#    plt.show()
#    plt.close()

def load_data():
    masks = np.load(OUTPUT_DIR + 'all_masks.npz')
    
#    X, _y = chooseAll(CSVFILE, 0, include_ID=True)

    x_test = [[],[],[],[],[]]
#    x_sep = [[],[],[],[],[]]
    
    filenames = glob.glob(PIXDIR + 'psc*.[3].fits')
    for full_filename in filenames:
        print(full_filename)
        dotSplit = full_filename.split('.')
        idNumString = dotSplit[-3].split('c')[-1] #used for getting hectospec data
        idNum  = int(idNumString)
        all_colors = []
        if len(glob.glob(full_filename[:-6] + '[3-6]' + full_filename[-5:])) !=4:
            print("missing filters for %s" % idNum)
            continue 
        for i in range(3,7):
            filename = full_filename[:-6] + str(i) + full_filename[-5:]   
            with fits.open(filename) as image_file:
                image_data = image_file[0].data
                exposure_time = float(image_file[0].header['EXPTIME'])
                image_data = image_data.byteswap(True).newbyteorder()
                # subtracting out background
                image_data = image_data - sep.Background(image_data) 
                image_data = image_data / exposure_time
        
            
            # find pixel coordinates of event
            event = db.where(db['eventID'] == idNum).dropna()
            eventRa = event['ra'].values[0] #values gives np arrays
            eventDec = event['dec'].values[0] #'hh:mm:ss.sss'
            
            # converting to degrees
            eventCoords = SkyCoord(eventRa, eventDec, unit=(u.hourangle, u.deg))
            eventRa = eventCoords.ra.deg
            eventDec = eventCoords.dec.deg
            #TODO: ircs or fk5?
            
            # get event pixel coords
            w = WCS(filename)
            eventX, eventY = w.all_world2pix(eventRa, eventDec, 1)
            #val = image_data[int(eventX), int(eventY)]
            paddedX = None
            paddedY = None
            if size(0, image_data) != 240:
                if eventX < 120:
                    paddedX = 'above'
                    image_data = pad(0, 'above', image_data)
                else:
                    paddedX = 'below'
                    image_data = pad(0, 'below', image_data)
            if size(1, image_data) != 240:
                if eventY < 120:
                    paddedY = 'above'
                    image_data = pad(1, 'above', image_data)
                else:
                    paddedY = 'below'
                    image_data = pad(1, 'below', image_data)
                    
            if size(0, image_data) != 240 or size(1, image_data) != 240:
                print('width, height: %s, %s' % (size(0, image_data), size(1, image_data)))
                raise Exception("bad dimensions")
            all_colors.append(image_data)    
            #plot(image_data, eventX, eventY)
        shp = image_data.shape
        redshift = zdict[idNumString]
        all_colors.append(np.full(shp, redshift))

        try:
            mask = masks[idNumString]
        except KeyError:
            print("No mask found for %s in all_masks.npz" % idNumString)
            continue
        if paddedX:
            mask = pad(0, paddedX, mask)
        if paddedY:
            mask = pad(1, paddedY, mask)
        if mask.shape != all_colors[0].shape:
            raise
        all_colors.append(mask)

        all_colors = np.nan_to_num(np.array(all_colors))
        typ = intDict[typeDict[idNumString]]
        x_test[typ].append(all_colors)

#        #lookup IDnum in SEP properties table from csv
#        rows, columns = np.asarray(X==idNum).nonzero()
#        added = False
#
#
#        for i in range(len(rows)):
#            if columns[i] == 0: # match is in the ID columnd
#                props = X[rows[i]]
#                x_sep[typ].append(props)
#                added = True
#                break
#        if not added:
#            print("no properties found for " + idNumString)

    for i in range(len(x_test)):
        x_test[i] = np.array(x_test[i])
        x_test[i] = np.transpose(x_test[i], [0,2,3,1])
        np.save(OUTPUT_DIR + "x_all2_%s" % i, x_test[i])
#        np.save(OUTPUT_DIR + "x_sep2_%s" % i, x_sep[i])

def main():
    load_data()
    
if __name__ == "__main__":
     main()
